"""Render the formula schematic and reference PDF from the package's Rd text."""
from pathlib import Path
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.font_manager import FontProperties

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs"
font_path = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
font = FontProperties(fname=font_path)
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"],
                     "mathtext.fontset": "stix", "svg.fonttype": "none"})
fig, ax = plt.subplots(figsize=(14, 10))
fig.subplots_adjust(left=.02, right=.98, top=.98, bottom=.02)
ax.set(xlim=(0, 14), ylim=(0, 10)); ax.axis("off")
ink, blue, teal, rust = "#222222", "#34618b", "#267d79", "#a25155"

def box(x, y, w, h, color, title, lines):
    ax.add_patch(Rectangle((x, y), w, h, facecolor="white", edgecolor=color, lw=1.3))
    ax.text(x+w/2, y+h-.22, title, ha="center", va="top", color=color,
            fontsize=14, fontweight="bold")
    for dy, value, size in lines:
        ax.text(x+w/2, y+dy, value, ha="center", va="center", fontsize=size, color=ink)

def arrow(x, y, xx, yy, color=blue):
    ax.add_patch(FancyArrowPatch((x, y), (xx, yy), arrowstyle="-|>",
                                mutation_scale=14, color=color, lw=1.2))

ax.text(7, 9.86, "Weighted Beta-mixture updates for BOLD extensions",
        ha="center", va="top", fontsize=21, fontweight="bold")
ax.text(7, 9.38, "Formula schematic of the implemented sequential scoring rule",
        ha="center", fontsize=12, color="#555555")
box(2.3, 7.65, 9.4, 1.35, blue, "Parent component after resolved outcomes", [
    (.73, r'$p_j\mid h\sim\mathrm{Beta}(A_h,B_h),\qquad w_h\geq0,\quad\sum_h w_h=1$', 18),
    (.28, r'$A_h=a_j+y_j+k_h,\qquad B_h=b_j+n_j^{\mathrm{obs}}-y_j+m_h-k_h$', 16)])
arrow(5, 7.65, 3.5, 7.12); arrow(9, 7.65, 10.5, 7.12, teal)
box(.25, 4.65, 6.55, 2.48, rust, "Hypothetical untreated DLT: z = 1", [
    (1.77, r'$\mathrm{Beta}(A_h+1,B_h),\qquad\theta_{h1}=\frac{A_h+1}{A_h+B_h+1}$', 18),
    (.97, r'$s_{h1}=\frac{\ell_1 q\,\theta_{h1}}{\ell_1 q\,\theta_{h1}+\ell_0(1-\theta_{h1})}$', 22),
    (.28, r'$r_{h1}=w_hs_{h1}$', 18)])
box(7.2, 4.65, 6.55, 2.48, teal, "Hypothetical untreated non-DLT: z = 0", [
    (1.77, r'$\mathrm{Beta}(A_h,B_h+1),\qquad\theta_{h0}=\frac{A_h}{A_h+B_h+1}$', 18),
    (.97, r'$s_{h0}=\frac{\ell_0(1-\theta_{h0})}{\ell_1 q\,\theta_{h0}+\ell_0(1-\theta_{h0})}$', 22),
    (.28, r'$r_{h0}=w_hs_{h0}$', 18)])
arrow(3.5, 4.65, 5.5, 4.12, rust); arrow(10.5, 4.65, 8.5, 4.12, teal)
box(1.1, 2.65, 11.8, 1.5, blue, "Normalize globally, then merge identical Beta components", [
    (.85, r'$\widetilde w_{hz}=\frac{r_{hz}}{\sum_{h^{\prime}}\sum_{z^{\prime}=0}^{1}r_{h^{\prime}z^{\prime}}},\qquad f_{\mathrm{new}}(p_j)=\sum_h\sum_{z=0}^{1}\widetilde w_{hz}\,\mathrm{Beta}(p_j;A_h+z,B_h+1-z)$', 18),
    (.28, "Repeat for each latent patient; equal-count branches share one component.", 12)])
ax.text(.4, 2.2, r'Prediction: $L=+\Rightarrow(\ell_1,\ell_0)=(s,1-c)$; '
        r'$L=-\Rightarrow(1-s,c)$; no prediction $\Rightarrow(1,1)$.', fontsize=14)
ax.text(.4, 1.78, r'Pending untreated outcome: $q=1-t/W$. Intervention without observed DLT: $q=1$.', fontsize=14)
ax.text(.4, 1.36, r'Observed intervention failure: ordinary DLT update $(A_h,B_h)\mapsto(A_h+1,B_h)$; no new branch.', fontsize=14)
ax.text(.4, .86, r'$s$: sensitivity; $c$: specificity; $W$: DLT window. Prior: '
        r'$a_j=\mu_j\nu_j,\ b_j=(1-\mu_j)\nu_j$.', fontsize=13)
for ext in ("png", "svg", "pdf"):
    fig.savefig(DOC / f"beta-mixture-formulas.{ext}", dpi=240, facecolor="white")
plt.close(fig)

# Generate the conventional R manual, including package metadata and index.
subprocess.run(["R", "CMD", "Rd2pdf", "--no-preview", "--force",
                "--output=docs/BOLDx-reference-manual.pdf", "."], cwd=ROOT, check=True)
