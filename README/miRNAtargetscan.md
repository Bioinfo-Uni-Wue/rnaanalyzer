
# Recommended IntaRNA Parameters for miRNA Target Prediction

This document outlines optimal parameters to use when applying **IntaRNA** for scanning a miRNA database against user-provided sequences in a web application environment.

---

## 🔧 Key IntaRNA Parameters

### 1. Seed Region Constraints

- **`--seedBP`**: Minimum number of base pairs in the seed region.
  - **Recommended**: `7`
- **`--seedMaxUP`**: Maximum number of unpaired bases allowed in the seed region.
  - **Recommended**: `1`

> These settings enforce a canonical 7-mer seed match, crucial for miRNA-mRNA recognition.

---

### 2. Accessibility Constraints

- **`--accW`**: Window size for accessibility computation.
  - **Recommended**: `150`
- **`--accL`**: Length of target region used for accessibility calculation.
  - **Recommended**: `100`

> Accessibility settings ensure local RNA structure is considered, reflecting true biological availability for interaction.

---

### 3. Interaction Lengths

- **`--qIntLenMax`**: Max interaction length on the query (miRNA).
  - **Recommended**: `25`
- **`--tIntLenMax`**: Max interaction length on the target (e.g., mRNA).
  - **Recommended**: `25`

> Keeps predictions biologically realistic by limiting interaction spans.

---

### 4. Energy Thresholds

- **`--outFilterEnergy`**: Minimum interaction energy to consider.
  - **Recommended**: `-15`

> Filters out weak, non-functional interactions.

---

### 5. Output Format

- **`--outMode`**: Output format.
  - **Recommended**: `C` (CSV format)

> Useful for structured, parsable results in a web application.

---

## 📌 Example Usage

```bash
IntaRNA \
  --query miRNA_database.fa \
  --target user_sequence.fa \
  --seedBP 7 \
  --seedMaxUP 1 \
  --accW 150 \
  --accL 100 \
  --qIntLenMax 25 \
  --tIntLenMax 25 \
  --outFilterEnergy -15 \
  --outMode C
```

---

## 📚 References

1. Mann, M., Wright, P. R., & Backofen, R. (2017). IntaRNA 2.0: enhanced and customizable prediction of RNA–RNA interactions. *Nucleic Acids Research*, 45(W1), W435–W439.
2. ViennaRNA Package Documentation: https://www.tbi.univie.ac.at/RNA/
3. IntaRNA GitHub Repository: https://github.com/BackofenLab/IntaRNA

---

For questions or integration guidance, please contact your system administrator or refer to the official IntaRNA documentation.
