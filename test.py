from transformers import AutoTokenizer, T5ForConditionalGeneration

# -------- Model (CPU-only) --------
MODEL_ID = "QizhiPei/biot5-base"  # you can try "QizhiPei/biot5-plus-base" if drift persists
tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, model_max_length=512)
model = T5ForConditionalGeneration.from_pretrained(MODEL_ID)

# -------- Prompt components --------
FEW_SHOT = """\
INSTRUCTION:
You are a biomedical assistant. Summarize RNA analysis results in exactly three titled sections.

FORMAT:
1) Interpretive summary
2) What this suggests
3) Recommended follow-ups

EVIDENCE:
Highly structured RNA. IRE motif at ~nt 34. No credible ORF. Several RNA-binding proteins detected.

OUTPUT:
1) Interpretive summary
The RNA forms a stable hairpin typical of cis-regulatory elements. The iron-responsive element near ~34 nt implicates iron homeostasis via IRP binding. Lack of coding features and the presence of multiple RBPs point to post-transcriptional control rather than translation.
2) What this suggests
Noncoding regulatory RNA likely modulating stability/translation in response to iron.
3) Recommended follow-ups
Validate IRP binding (EMSA), test iron-dependent reporter activity, assess conservation of the IRE stem-loop, and perform targeted RBP pulldown for the detected binders.
"""

INSTRUCTION = """\
INSTRUCTION:
You are a biomedical assistant. Summarize RNA analysis results in exactly three titled sections.

FORMAT:
1) Interpretive summary
2) What this suggests
3) Recommended follow-ups
"""

def build_prompt(evidence: str) -> str:
    # T5 benefits from a task prefix; keep it concise and deterministic.
    return "summarize: " + FEW_SHOT + "\n" + INSTRUCTION + "\nEVIDENCE:\n" + evidence.strip() + "\n[END]"

# -------- Decoding with strong anti-looping & forced headings --------
def generate_summary(evidence: str, max_new_tokens: int = 220) -> str:
    prompt = build_prompt(evidence)
    inputs = tokenizer(prompt, return_tensors="pt", truncation=True)

    # Ban common HTML-ish junk & separators that caused your logs.
    junk = [
        "<p>", "</p>", "<span>", "</span>", "<div>", "</div>",
        "<br>", "</br>", "<", ">", "</", "<P>", "</P>"
    ]
    bad_words_ids = []
    for s in junk:
        ids = tokenizer.encode(s, add_special_tokens=False)
        if ids:
            bad_words_ids.append(ids)

    # Force the three section titles to appear (helps anchor generation).
    force_words = [
        tokenizer.encode("1) Interpretive summary", add_special_tokens=False),
        tokenizer.encode("2) What this suggests", add_special_tokens=False),
        tokenizer.encode("3) Recommended follow-ups", add_special_tokens=False),
    ]

    out = model.generate(
        **inputs,
        max_new_tokens=max_new_tokens,
        num_beams=6,
        length_penalty=0.9,
        no_repeat_ngram_size=4,
        repetition_penalty=1.25,
        early_stopping=True,
        eos_token_id=tokenizer.eos_token_id,
        bad_words_ids=bad_words_ids if bad_words_ids else None,
        force_words_ids=force_words
    )
    return tokenizer.decode(out[0], skip_special_tokens=True)

# -------- Example usage --------
example_evidence = """Highly structured RNA present; likely regulatory function.
IRE motif at position 34; potential IRE function.
noncoding RNA.
3 RNA-binding proteins detected.
"""

if __name__ == "__main__":
    print(generate_summary(example_evidence))

    print("\n--- Sample 2 ---")
    print(generate_summary("Low structure; strong 420aa ORF; Ribo-seq footprints; no known motifs."))

