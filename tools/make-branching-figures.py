"""Two formula-only branching schematics matching the BOLDx observation rules."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle

DOC = Path(__file__).resolve().parents[1] / "docs"
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"],
                     "mathtext.fontset": "stix", "svg.fonttype": "none"})
INK = "#222a31"
BLUE = "#355f89"
TEAL = "#147d83"
RED = "#b04d59"
GOLD = "#a57b28"
GREY = "#63727c"
PALE = "#f4f7f8"

def setup(title, subtitle):
    fig, ax = plt.subplots(figsize=(16, 14))
    fig.subplots_adjust(left=.02, right=.98, bottom=.015, top=.985)
    ax.set(xlim=(0,16), ylim=(0,14)); ax.axis("off")
    ax.text(8,13.75,title,fontsize=24,fontfamily="STIXGeneral",fontstyle="italic",
            fontweight="bold",color=INK,ha="center",va="top")
    ax.text(.35,13.22,subtitle,fontsize=13,color=GREY,va="top")
    return fig, ax

def panel(ax,x,y,w,h,title,color=BLUE,fill="white",lines=(),size=15):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.03,rounding_size=0.06",
                 lw=1.1,edgecolor=color,facecolor=fill,zorder=2))
    ax.text(x+w/2,y+h-.15,title,ha="center",va="top",fontsize=size,
            fontweight="bold",color=color,zorder=3)
    for dy,text,fs in lines:
        ax.text(x+w/2,y+dy,text,ha="center",va="center",fontsize=fs,color=INK,zorder=3)

def arrow(ax,start,end,color=GREY,curve=0,dashed=False):
    ax.add_patch(FancyArrowPatch(start,end,arrowstyle="-|>",mutation_scale=13,
                 connectionstyle=f"arc3,rad={curve}",lw=1.3,color=color,
                 linestyle="--" if dashed else "-",zorder=1))

def person(ax,x,y,color,label=None):
    ax.add_patch(Circle((x,y+.18),.075,facecolor=color,edgecolor="white",lw=.5,zorder=3))
    ax.add_patch(FancyBboxPatch((x-.085,y-.08),.17,.18,
        boxstyle="round,pad=0.01,rounding_size=0.04",facecolor=color,edgecolor=color,zorder=3))
    ax.plot([x-.04,x-.04],[y-.08,y-.19],color=color,lw=2.4)
    ax.plot([x+.04,x+.04],[y-.08,y-.19],color=color,lw=2.4)
    if label:
        ax.text(x+.17,y+.08,label,fontsize=11,color=color,va="center")

def clock(ax,x,y,color=BLUE):
    ax.add_patch(Circle((x,y),.15,facecolor="white",edgecolor=color,lw=1.3))
    ax.plot([x,x,x+.08],[y+.09,y,y-.03],color=color,lw=1.3)

def common_formulas(ax,top=4.25):
    ax.plot([.35,15.65],[top,top],color="#d3dde2",lw=1)
    ax.text(.35,top-.18,"Shared weighted Beta-mixture update",fontsize=16,
            fontweight="bold",color=BLUE,va="top")
    ax.text(.35,top-.65,
        r'$A_h=a_j+y_j+k_h,\quad B_h=b_j+n_j^{\mathrm{obs}}-y_j+m_h-k_h;'
        r'\qquad p_j\mid h\sim\mathrm{Beta}(A_h,B_h),\quad w_h\geq0,\quad\sum_h w_h=1$',fontsize=16,color=INK)
    panel(ax,.35,top-2.35,7.45,1.35,r'Untreated DLT child: $z=1$',RED,"#fcf5f6",[
        (.86,r'$\mathrm{Beta}(A_h+1,B_h),\quad\theta_{h1}=\frac{A_h+1}{A_h+B_h+1}$',16),
        (.35,r'$s_{h1}=\frac{\ell_1q\,\theta_{h1}}{\ell_1q\,\theta_{h1}+\ell_0(1-\theta_{h1})},\quad r_{h1}=w_hs_{h1}$',17)],size=13)
    panel(ax,8.15,top-2.35,7.5,1.35,r'Untreated non-DLT child: $z=0$',TEAL,"#f1f8f7",[
        (.86,r'$\mathrm{Beta}(A_h,B_h+1),\quad\theta_{h0}=\frac{A_h}{A_h+B_h+1}$',16),
        (.35,r'$s_{h0}=\frac{\ell_0(1-\theta_{h0})}{\ell_1q\,\theta_{h0}+\ell_0(1-\theta_{h0})},\quad r_{h0}=w_hs_{h0}$',17)],size=13)
    ax.text(8,top-2.95,
        r'$\widetilde w_{hz}=\frac{r_{hz}}{\sum_{h^{\prime}}\sum_{z^{\prime}=0}^{1}r_{h^{\prime}z^{\prime}}},'
        r'\qquad f_{\mathrm{new}}(p_j)=\sum_h\sum_{z=0}^{1}\widetilde w_{hz}\,'
        r'\mathrm{Beta}(p_j;A_h+z,B_h+1-z)$',ha="center",fontsize=17,color=INK)
    ax.text(.35,top-3.6,
        r'$L=+:(\ell_1,\ell_0)=(s,1-c);\quad L=-:(1-s,c);\quad '
        r'\mathrm{no\ label}:(1,1).\qquad a_j=\mu_j\nu_j,\quad b_j=(1-\mu_j)\nu_j.$',
        fontsize=14,color=INK)
    ax.text(.35,top-4.0,
        r'$s$: sensitivity; $c$: specificity; $W$: DLT window; $w_h$: parent weight; '
        r'$m_h$: prior latent patients; $k_h$: their hypothesized DLT count.',fontsize=11,color=GREY)

def save(fig,name):
    for ext in ("png","pdf","svg"):
        fig.savefig(DOC/f"{name}.{ext}",dpi=260,facecolor="white")
    plt.close(fig)

# Variant 1: two positive predictions followed by preventive intervention.
fig, ax = setup("BOLD with prediction-guided intervention",
    "Two positive predictions at one dose; no observed intervention failure. Untreated outcomes remain latent.")
for x, label, col in [(1.8,"DOSE-SPECIFIC PRIOR",BLUE),(6.05,"AFTER PATIENT 1",GOLD),(12.1,"AFTER PATIENT 2",GOLD)]:
    ax.text(x,12.65,label,fontsize=12,ha="center",color=col,fontweight="bold")
person(ax,5.2,12.25,GOLD,r'$L_1=+$'); person(ax,11.25,12.25,GOLD,r'$L_2=+$')
panel(ax,.4,8.05,2.85,1.5,"Start at dose j",BLUE,PALE,[
    (.81,r'$\mathrm{Beta}(a_j,b_j)$',19),(.32,r'$w_{\varnothing}=1$',18)])
panel(ax,4.25,9.85,3.6,1.45,r'$z_1=1$',RED,"#fcf5f6",[
    (.84,r'$\mathrm{Beta}(a_j+1,b_j)$',18),(.29,r'$w_1^{(1)}=\frac{r_1^{(1)}}{r_0^{(1)}+r_1^{(1)}}$',18)])
panel(ax,4.25,6.25,3.6,1.45,r'$z_1=0$',TEAL,"#f1f8f7",[
    (.84,r'$\mathrm{Beta}(a_j,b_j+1)$',18),(.29,r'$w_0^{(1)}=\frac{r_0^{(1)}}{r_0^{(1)}+r_1^{(1)}}$',18)])
arrow(ax,(3.25,8.8),(4.25,10.55),RED,-.1)
arrow(ax,(3.25,8.8),(4.25,6.98),TEAL,.1)
leaves=[("11",(1,1),11.0,RED,r'\mathrm{Beta}(a_j+2,b_j)'),
        ("10",(1,0),9.25,TEAL,r'\mathrm{Beta}(a_j+1,b_j+1)'),
        ("01",(0,1),7.5,RED,r'\mathrm{Beta}(a_j+1,b_j+1)'),
        ("00",(0,0),5.75,TEAL,r'\mathrm{Beta}(a_j,b_j+2)')]
for path,(u,v),y,col,beta in leaves:
    panel(ax,10.0,y-.4,5.55,1.43,fr'Path {path}: $(z_1,z_2)=({u},{v})$',col,
        "#fcf5f6" if v else "#f1f8f7",[
        (.82,f'${beta}$',18),
        (.28,fr'$v_{{{path}}}=w_{{{u}}}^{{(1)}}s_{{{path}}}^{{(2)}},\quad '
              fr'\widehat w_{{{path}}}=v_{{{path}}}/\sum_{{u,v}}v_{{uv}}$',16)],size=13)
    arrow(ax,(7.85,10.55 if u else 6.98),(10,y+.3),col,
        .06 if v else -.06)
ax.text(.55,7.1,"Each positive prediction splits\nevery current component.",
        fontsize=12,color=GREY,ha="left",va="center")
clock(ax,1.0,5.9,GOLD)
ax.text(1.3,5.9,r'$q=1$: no-DLT follow-up does not reveal untreated status.',
        fontsize=13,color=INK,va="center")
ax.text(.5,5.23,r'$f_2(p_j)=\widehat w_{11}\mathrm{Beta}(p_j;a_j+2,b_j)'
        r'+(\widehat w_{10}+\widehat w_{01})\mathrm{Beta}(p_j;a_j+1,b_j+1)'
        r'+\widehat w_{00}\mathrm{Beta}(p_j;a_j,b_j+2)$',fontsize=17,color=INK)
ax.text(.5,4.65,"Normalize all four paths together; merge the two equal-parameter components. "
        "Observed failure instead enters as a known DLT.",fontsize=12,color=GREY)
common_formulas(ax)
save(fig,"prediction-branching-schematic")

# Variant 2: temporal resolution for observed untreated outcomes, with the
# intervention counterfactual exception shown as a separate parallel lane.
fig, ax = setup("TITE-BOLD with prediction-guided intervention",
    "Pending untreated outcomes are replaced by observations; prevented-toxicity counterfactuals remain latent.")
ax.plot([1,15.1],[12.35,12.35],color="#d5dfe4",lw=2)
for x,lab in [(1.1,"ENROLLMENT"),(6.2,"PARTIAL FOLLOW-UP"),(12.55,"OUTCOME RESOLUTION")]:
    clock(ax,x,12.35,BLUE)
    ax.text(x,11.93,lab,ha="center",fontsize=11,color=BLUE,fontweight="bold")

ax.text(.4,11.35,"A  No preventive intervention",fontsize=16,color=BLUE,fontweight="bold")
person(ax,.65,10.65,BLUE)
panel(ax,1.02,9.85,3.0,1.22,"Prediction + pending outcome",BLUE,PALE,[
    (.61,r'$p_j\mid h\sim\mathrm{Beta}(A_h,B_h)$',15),
    (.22,r'$q(t)=1-t/W$',17)],size=11)
panel(ax,4.8,10.55,3.65,.9,r'$z_i=1$',RED,"#fcf5f6",[
    (.34,r'$\mathrm{Beta}(A_h+1,B_h),\ r_{h1}(t)$',15)],size=12)
panel(ax,4.8,9.3,3.65,.9,r'$z_i=0$',TEAL,"#f1f8f7",[
    (.34,r'$\mathrm{Beta}(A_h,B_h+1),\ r_{h0}(t)$',15)],size=12)
arrow(ax,(4.05,10.47),(4.8,11),RED,-.15)
arrow(ax,(4.05,10.47),(4.8,9.75),TEAL,.15)
panel(ax,10.0,10.55,5.55,.9,r'Observed DLT: $Y_i=1$',RED,"#fcf5f6",[
    (.33,r'$(A_h,B_h)\longmapsto(A_h+1,B_h)$',17)],size=13)
panel(ax,10.0,9.3,5.55,.9,r'No DLT through $W$: $Y_i=0$',TEAL,"#f1f8f7",[
    (.33,r'$(A_h,B_h)\longmapsto(A_h,B_h+1)$',17)],size=13)
arrow(ax,(8.48,11),(10,11),RED)
arrow(ax,(8.48,9.75),(10,9.75),TEAL)
ax.text(10.0,8.95,"Replace this patient's pending label/scores; rebuild from observed history.",
        fontsize=11,color=BLUE)

ax.plot([.4,15.6],[8.65,8.65],color="#dbe3e6",lw=1)
ax.text(.4,8.3,"B  Preventive intervention after a positive prediction",fontsize=16,color=GOLD,fontweight="bold")
person(ax,.65,7.55,GOLD)
panel(ax,1.02,6.5,3.0,1.35,r'$L_i=+$; intervene',GOLD,"#fbf8ef",[
    (.75,r'$q=1$',19),(.25,"Untreated status is unknown",11)],size=13)
panel(ax,4.8,7.3,3.65,.9,r'$z_i=1$',RED,"#fcf5f6",[
    (.34,r'$\mathrm{Beta}(A_h+1,B_h),\ r_{h1}$',15)],size=12)
panel(ax,4.8,6.05,3.65,.9,r'$z_i=0$',TEAL,"#f1f8f7",[
    (.34,r'$\mathrm{Beta}(A_h,B_h+1),\ r_{h0}$',15)],size=12)
arrow(ax,(4.05,7.17),(4.8,7.75),RED,-.15)
arrow(ax,(4.05,7.17),(4.8,6.5),TEAL,.15)
panel(ax,10,7.3,5.55,.9,"Observed intervention failure",RED,"#fcf5f6",[
    (.34,r'$Y_i=1:\quad(A_h,B_h)\mapsto(A_h+1,B_h)$',17)],size=13)
panel(ax,10,5.85,5.55,1.1,"No observed DLT after intervention",GOLD,"#fbf8ef",[
    (.6,r'$Y_i=0\ \not\Rightarrow\ z_i=0$',18),
    (.21,r'Both untreated-outcome branches remain; $q=1$.',12)],size=12)
arrow(ax,(8.48,7.75),(10,7.75),RED)
arrow(ax,(8.48,7.75),(10,6.4),GOLD,.12,dashed=True)
arrow(ax,(8.48,6.5),(10,6.4),GOLD)
ax.text(.5,5.35,r'Observed outcomes replace predictions only when untreated status becomes known. '
        r'After prevention, the unobserved counterfactual does not resolve.',fontsize=12,color=INK)
ax.text(.5,4.75,"Current combined protocol: intervened patients finish follow-up before the next dose decision; "
        "other outcomes may remain pending.",fontsize=12,color=GREY)
common_formulas(ax)
save(fig,"tite-prediction-branching-schematic")
