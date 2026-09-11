"""
Step 4/5: PIK3CA sliding-window handling + ESM2 masked-marginal scoring.

For proteins within the 1022-residue limit: score directly on the full
canonical sequence.

For PIK3CA (1068 aa, exceeds limit): extract a window of up to 1023 residues
centered on the mutated position, shifted inward near either terminus so the
window never runs off the sequence. The mutation's position is re-indexed
relative to the window before masking/scoring. context_window_applied=True
is recorded on every PIK3CA output row.

Includes: model-load verification, a determinism canary (score a fixed
variant twice, confirm identical output), and per-variant scoring.
"""
import json
import os
import sys
import time
import torch
from transformers import AutoTokenizer, EsmForMaskedLM
from pipeline_utils import get_window as _get_window, ESM2_CONTEXT_LIMIT

CHECKPOINT_PATH = "data/scored/_checkpoint.json"
MAX_SECONDS_PER_RUN = 150  # stay safely under the 300s tool limit, allowing for slow loads

MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
EXPECTED_PARAMS = 651_000_000  # approx, from feasibility check
CONTEXT_LIMIT = ESM2_CONTEXT_LIMIT
WINDOW_SIZE = 1023  # per locked decision: up to 1023 residues centered on the mutation


def get_window(seq, pos, window_size=WINDOW_SIZE):
    """Wraps the shared get_window() (src/pipeline_utils.py) to additionally
    report whether a window was actually applied, for output flagging."""
    windowed_seq, new_pos = _get_window(seq, pos, window_size)
    window_applied = len(seq) > window_size
    return windowed_seq, new_pos, window_applied


def load_model():
    print(f"Loading {MODEL_NAME}...")
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForMaskedLM.from_pretrained(MODEL_NAME)
    model.eval()
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  Loaded in {time.time()-t0:.1f}s | {n_params/1e6:.1f}M params | "
          f"max_position_embeddings={model.config.max_position_embeddings}")
    # model-load sanity check
    assert abs(n_params - EXPECTED_PARAMS) / EXPECTED_PARAMS < 0.01, \
        f"Loaded model param count ({n_params}) doesn't match expected checkpoint!"
    assert model.config.max_position_embeddings == 1026, \
        "Context limit assumption violated - re-verify before scoring!"
    return tokenizer, model


def score_variant(seq, pos, wt, mt, tokenizer, model):
    """masked-marginal score = logP(mt) - logP(wt) at the masked position."""
    inputs = tokenizer(seq, return_tensors="pt")
    tokens = inputs["input_ids"].clone()
    token_pos = pos  # CLS occupies index 0, so residue i (1-indexed) -> token index i
    assert tokens.shape[1] > token_pos, "position falls outside tokenized sequence!"

    original_token = tokenizer.convert_ids_to_tokens([tokens[0, token_pos].item()])[0]
    assert original_token == wt, (
        f"Tokenized sequence residue at masked position ({original_token}) "
        f"doesn't match expected WT ({wt}) - indexing bug!"
    )

    tokens[0, token_pos] = tokenizer.mask_token_id
    with torch.no_grad():
        logits = model(input_ids=tokens, attention_mask=inputs["attention_mask"]).logits
    log_probs = torch.log_softmax(logits[0, token_pos], dim=-1)
    wt_id = tokenizer.convert_tokens_to_ids(wt)
    mt_id = tokenizer.convert_tokens_to_ids(mt)
    return (log_probs[mt_id] - log_probs[wt_id]).item()


def main():
    with open("data/processed/unique_variants_to_score.json") as f:
        variants = json.load(f)
    with open("data/raw/sequences.json") as f:
        sequences = json.load(f)

    # ---- resume from checkpoint if one exists ----
    if os.path.exists(CHECKPOINT_PATH):
        with open(CHECKPOINT_PATH) as f:
            checkpoint = json.load(f)
        results = checkpoint["results"]
        start_idx = checkpoint["next_idx"]
        print(f"Resuming from checkpoint: {start_idx}/{len(variants)} already scored.")
    else:
        results = []
        start_idx = 0

    if start_idx >= len(variants):
        print("All variants already scored. Finalizing.")
        finalize(results)
        return

    tokenizer, model = load_model()

    # ---- determinism canary: score one fixed variant twice (only on first run) ----
    if start_idx == 0:
        canary = next(v for v in variants if v["gene"] == "TP53")
        seq = sequences["TP53"]["sequence"]
        s1 = score_variant(seq, canary["protein_position"], canary["wt_residue"],
                            canary["mutant_residue"], tokenizer, model)
        s2 = score_variant(seq, canary["protein_position"], canary["wt_residue"],
                            canary["mutant_residue"], tokenizer, model)
        assert s1 == s2, f"DETERMINISM CHECK FAILED: {s1} != {s2}"
        print(f"Determinism canary passed (TP53 {canary['variant_notation']}: {s1:.4f} both runs)\n")

    print(f"Scoring variants {start_idx}..{len(variants)-1} (batch time budget {MAX_SECONDS_PER_RUN}s)...")
    t_batch_start = time.time()
    per_gene_window_count = {}
    i = start_idx

    for i in range(start_idx, len(variants)):
        if time.time() - t_batch_start > MAX_SECONDS_PER_RUN:
            print(f"\nTime budget reached at index {i}. Checkpointing and stopping this batch.")
            break
        v = variants[i]
        gene = v["gene"]
        seq_info = sequences[gene]
        full_seq = seq_info["sequence"]
        pos = v["protein_position"]
        wt = v["wt_residue"]
        mt = v["mutant_residue"]

        # sanity re-check against full sequence before windowing
        assert full_seq[pos - 1] == wt, f"Pre-scoring WT check failed for {gene} {v['variant_notation']}!"

        windowed_seq, window_pos, window_applied = get_window(full_seq, pos)
        if window_applied:
            assert windowed_seq[window_pos - 1] == wt, (
                f"Windowing re-indexing bug for {gene} {v['variant_notation']}!"
            )
            per_gene_window_count[gene] = per_gene_window_count.get(gene, 0) + 1

        score = score_variant(windowed_seq, window_pos, wt, mt, tokenizer, model)

        results.append({
            **v,
            "esm2_score": score,
            "esm2_model": MODEL_NAME,
            "context_window_applied": window_applied,
            "scored_sequence_length": len(windowed_seq),
        })

        if (i + 1 - start_idx) % 10 == 0:
            elapsed = time.time() - t_batch_start
            rate = elapsed / (i + 1 - start_idx)
            print(f"  {i+1}/{len(variants)} scored total ({rate:.2f}s/variant this batch)")
            # incremental checkpoint - never lose more than ~10 variants to an external kill
            with open(CHECKPOINT_PATH, "w") as f:
                json.dump({"results": results, "next_idx": i + 1}, f)

    next_idx = i + 1 if (time.time() - t_batch_start <= MAX_SECONDS_PER_RUN) else i

    if per_gene_window_count:
        print(f"Context-window applied this batch: {per_gene_window_count}")

    # checkpoint progress
    with open(CHECKPOINT_PATH, "w") as f:
        json.dump({"results": results, "next_idx": next_idx}, f)
    print(f"\nCheckpoint saved: {next_idx}/{len(variants)} variants scored so far.")

    if next_idx >= len(variants):
        print("All variants scored. Finalizing.")
        finalize(results)
    else:
        print(f"Run this script again to continue from index {next_idx}.")


def finalize(results):
    scores = [r["esm2_score"] for r in results]
    print(f"\nFinal scored set: {len(results)} unique variants")
    print(f"Score range: min={min(scores):.2f}, max={max(scores):.2f}")
    out_of_range = [r for r in results if not (-20 <= r["esm2_score"] <= 5)]
    if out_of_range:
        print(f"  NOTE: {len(out_of_range)} scores outside expected [-20,5] plausibility band - flagged for review:")
        for r in out_of_range:
            print(f"    {r['gene']} {r['variant_notation']}: {r['esm2_score']:.2f}")

    with open("data/scored/scored_variants.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved data/scored/scored_variants.json")
    if os.path.exists(CHECKPOINT_PATH):
        os.remove(CHECKPOINT_PATH)


if __name__ == "__main__":
    main()
