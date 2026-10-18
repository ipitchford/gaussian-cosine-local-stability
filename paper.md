# Local stability of cosine configurations\
for Gaussian minima

Anonymous

18 October 2026

## Abstract

For a centred Gaussian vector with unit coordinate variances, how small can the minimum absolute coordinate be? The equally spaced planar cosine configuration is the proposed global minimizer, but the higher-rank comparison remains unresolved. We prove that for every number of coordinates $`n\ge3`$ and every moment order $`p>1`$, this configuration is a strict local minimum in the full correlation domain, including perturbations that increase covariance rank. An explicit second-variation formula separates angular perturbations from normal perturbations. The normal form is a positive multiple of a shifted signed-cycle adjacency matrix; its kernel consists exactly of ambient rotations. We justify continuous second derivatives across moving active interfaces, including the integrable zero-coordinate singularity for $`1<p<2`$, and obtain quadratic growth in squared distance from the orthogonal orbit. For four coordinates and the second moment we also give a sharp independent-coordinate subclass, an analytic nonglobal local minimum at independence, and two finite evaluation formulas. An elementary independent-noise bound excludes an explicit interior region from possible counterexamples. The intrinsic-volume identity used in one evaluator is a specialization of known work. The global four-coordinate inequality and the stronger stochastic conjecture are not proved.

This accessible companion retains source labels; the PDF is authoritative for theorem and equation numbering.

<div class="center">

Version 1.0.0DOI: [10.5281/zenodo.23453252](https://doi.org/10.5281/zenodo.23453252)

</div>

**Keywords:** Gaussian minima; cosine covariance; second variation; local stability; signed cycle; intrinsic volumes.

# Introduction and main result

Let $`C`$ be a real $`n\times n`$ correlation matrix and $`G\sim N(0,C)`$. We study
``` math
F_p(C)=\mathbb E\min_{1\le i\le n}|G_i|^p,\qquad p>0.
```
A correlation matrix is positive semidefinite with diagonal entries one. Thus $`C=UU^{\mathsf T}`$ for a matrix whose rows $`u_i`$ are unit vectors, and $`G`$ has the law of $`UZ`$ for a standard Gaussian $`Z`$. The objective is also the expected $`p`$th power of the distance from $`Z`$ to the union of the central hyperplanes $`u_i^\perp`$. Row sign changes and permutations do not change it.

Kunisky proposes a stochastic comparison with the cosine frame
``` math
v_i=(\cos(i\pi/n),\sin(i\pi/n)),\quad 0\le i<n,
 \qquad C_{\rm cos}=VV^{\mathsf T}.
```
Its perpendicular lines divide the circle modulo $`\pi`$ into equal gaps. The frame has rank two, so a comparison restricted to planar configurations does not settle the conjecture. A competing configuration could leave the plane. The main result below rules out every sufficiently small such improvement, for all $`n\ge3`$ and $`p>1`$, rather than just perturbations preserving rank.

The cosine law and its values are known . Earlier Gaussian-minimum conjectures and comparisons appear in . Our contribution is the explicit full local quadratic form and its strictness after removing rotations. The planar and independent-block arguments are supporting comparisons. The intrinsic-volume identity in Section <a href="#sec:intrinsic" data-reference-type="ref" data-reference="sec:intrinsic">8</a> is attributed to Long ; no new general identity is claimed there. Related convex hulls of entire great subspheres are studied in the preliminary manuscript . Those are different bodies from our finite signed-vertex hulls, and its conditional volume comparison does not prove the Gaussian conjecture.

Fix $`n\ge3`$, $`d\ge2`$, and $`p>1`$. Embed the rows of $`V`$ in $`\mathbb R^d`$ by adding zero coordinates. Put
``` math
\alpha=\frac{\pi}{n},\quad h=\frac{\alpha}{2},\quad
 m_q=\mathbb ER^q=2^{q/2}\Gamma(1+q/2)\quad(q>-2),
```
where $`R`$ is the length of a two-dimensional standard Gaussian. Define the symmetric matrix $`A_n`$ by positive edges $`(i,i+1)`$ for $`0\le i<n-1`$ and a negative edge $`(0,n-1)`$, with all other entries zero. The restriction $`n\ge3`$ avoids overlapping edges in this definition.

<div id="thm:main" class="theorem">

**Theorem 1** (All-order transverse stability). *Let $`a_i\in\mathbb R`$, $`b_i\in\mathbb R^{d-2}`$, $`a_n=a_0`$, and let $`B`$ have rows $`b_i`$. Along the unit-row path
``` math
\begin{equation}
\label{eq:path}
 u_i(t)=\left(\sqrt{1-t^2\|b_i\|^2}\,v(i\alpha+ta_i),\,tb_i\right),
 \qquad v(\theta)=(\cos\theta,\sin\theta),
\end{equation}
```
one has
``` math
\begin{align}
 F_p(U(t)U(t)^{\mathsf T})
 &= F_p(C_{\rm cos})+t^2 Q_{n,p}(a,B)+o(t^2),\label{eq:main}\\
 Q_{n,p}(a,B)
 &=H_{n,p}\sum_{i=0}^{n-1}(a_{i+1}-a_i)^2
 +K_{n,p}\operatorname{tr}\!\left[B^{\mathsf T}(2\cos\alpha I_n-A_n)B\right],\label{eq:Q}\\
 H_{n,p}&=\frac{p m_p\sin(h)^{p-1}\cos h}{4\pi},\qquad
 K_{n,p}=\frac{p m_{p-2}\sin(h)^{p-1}}{4\pi\cos h}.\label{eq:HK}
\end{align}
```
For fixed $`n,d,p`$, the remainder is uniform over bounded tangent directions. The kernel consists exactly of $`a_i`$ constant and $`B=VM`$, with $`M\in\mathbb R^{2\times(d-2)}`$.*

*Let
``` math
D(U)=\inf_{Q\in O(d)}\sum_i\|u_i-Q(v_i,0)\|^2,\qquad
 \lambda_{n,p}=\min\{4H_{n,p}\sin^2\alpha,\,
            2K_{n,p}(\cos\alpha-\cos3\alpha)\}>0.
```
In a sufficiently small neighbourhood of the orthogonal orbit,
``` math
\begin{equation}
\label{eq:growth}
 F_p(UU^{\mathsf T})-F_p(C_{\rm cos})
 \ge \tfrac12\lambda_{n,p}D(U).
\end{equation}
```
Consequently $`C_{\rm cos}`$ is a strict local minimum among all correlation matrices, including those of higher rank, modulo the stated discrete symmetries.*

</div>

The displayed $`Q`$ is the coefficient of $`t^2`$; the Hessian evaluated twice on the velocity is $`2Q`$. In dimension two the normal directions are absent, and $`\lambda_{n,p}`$ is still a valid lower bound, not necessarily the exact smallest positive coefficient. No explicit neighbourhood radius, nor radius uniform in $`n,d,p`$, is asserted.

For $`p=2`$, $`H_{n,2}=\sin\alpha/(2\pi)`$ and $`K_{n,2}=\tan(\alpha/2)/(2\pi)`$. In particular at $`n=4`$, with $`P=I_4-\tfrac12VV^{\mathsf T}`$,
``` math
\begin{equation}
\label{eq:four}
 Q_{4,2}(a,B)=\frac{\sqrt2}{4\pi}\sum_i(a_{i+1}-a_i)^2
       +\frac{2-\sqrt2}{\pi}\|PB\|_F^2.
\end{equation}
```
The main theorem is local. It does not prove that a remote frame cannot have a smaller moment. Sections <a href="#sec:block" data-reference-type="ref" data-reference="sec:block">5</a>–<a href="#sec:noise" data-reference-type="ref" data-reference="sec:noise">9</a> give complementary four-coordinate comparisons; Sections <a href="#sec:interface" data-reference-type="ref" data-reference="sec:interface">7</a>–<a href="#sec:intrinsic" data-reference-type="ref" data-reference="sec:intrinsic">8</a> retain explicit global proof obligations rather than claiming to solve them.

# The planar comparison

Order the unoriented lines perpendicular to planar unit vectors $`u_i`$. Their successive gaps are $`\delta_i\ge0`$, summing to $`\pi`$; coincident lines have zero gaps.

<div id="prop:planar" class="proposition">

**Proposition 2**. *Equally spaced planar lines minimize $`\mathbb P(\min_i|u_i\cdot Z|>z)`$ for each $`z>0`$, and hence every positive moment. Equality requires equal gaps. In particular
``` math
\mathbb E\min_i(u_i\cdot Z)^2
 =1-\frac1\pi\sum_i\sin\delta_i
 \ge1-\frac n\pi\sin(\pi/n).
```*

</div>

<div class="proof">

*Proof.* Write $`Z=R(\cos\phi,\sin\phi)`$ with independent uniform angle. On a gap of length $`\delta`$, the distance to the nearer zero line is $`\min(t,\delta-t)`$. Thus
``` math
F_p=\frac{2m_p}{\pi}\sum_i\int_0^{\delta_i/2}\sin^p t\,dt,\qquad
 \mathbb P(\min_i|u_i\cdot Z|>z)
 =\frac2\pi\sum_i\int_0^{\delta_i/2}
   e^{-z^2/(2\sin^2t)}\,dt.
```
For positive $`z`$, the integrand is strictly increasing on $`(0,\pi/2)`$, so each gap integral is strictly convex in $`\delta`$. Jensen’s inequality proves the survival comparison and equality statement. Integrating survival functions proves the moment comparison. At $`p=2`$, direct integration gives the displayed formula. These arguments include zero gaps. ◻

</div>

# Regularity across moving interfaces

This section supplies the analytic justification for varying the frame, rather than differentiating only inside fixed winning cells. Write $`f_i(U,x)=|u_i\cdot x|^p`$. A winning cell is a region where one $`f_i`$ is smallest. Signed interfaces have equations
``` math
H_{ij,s}(U)=\{x:(u_i-su_j)\cdot x=0\},\qquad s\in\{-1,1\}.
```
Only the portions where $`|u_k\cdot x|\ge|u_i\cdot x|`$ for every $`k`$ are active.

<div id="lem:C2" class="lemma">

**Lemma 3**. *For $`p>1`$, $`U\mapsto\mathbb E\min_i|u_i\cdot Z|^p`$ is $`C^2`$ locally on the product of unit spheres whenever its unoriented rows are distinct. Along a smooth path, its second-order coefficient is the sum of the ordinary bulk coefficients on the winning cells and, for every unordered signed active interface, the contribution
``` math
\begin{equation}
\label{eq:interface-second}
 -\frac{p}{2}\int_{H_{ij,s}}
 \frac{|u_i\cdot x|^{p-1}}{\|u_i-su_j\|}
 ((\dot u_i-s\dot u_j)\cdot x)^2
 \mathbf1_{\rm active}\varphi_d(x)\,dH(x).
\end{equation}
```*

</div>

<div class="proof">

*Proof.* Away from ties and zero coordinates differentiation is ordinary differentiation of the unique active branch. Across a two-cell interface use $`\min(f_i,f_j)=(f_i+f_j-|f_i-f_j|)/2`$. Its second distributional derivative contributes $`-\delta(f_i-f_j)(\dot f_i-\dot f_j)^2`$. At $`u_i\cdot x=s\,u_j\cdot x=z\ne0`$,
``` math
\|\nabla_x(f_i-f_j)\|=p|z|^{p-1}\|u_i-su_j\|,\quad
 \dot f_i-\dot f_j
 =p|z|^{p-1}\operatorname{sign}(z)(\dot u_i-s\dot u_j)\cdot x.
```
Coarea and division by two for the Taylor coefficient give <a href="#eq:interface-second" data-reference-type="eqref" data-reference="eq:interface-second">[eq:interface-second]</a>. Equivalently this calculation follows by integrating on the cells bounded by the linear signed walls, avoiding a delta expression at their intersections.

Here are the domination and null-set details needed to extend the calculation. Pair separations remain bounded away from zero in a small frame neighbourhood. Each interface therefore has a smoothly varying orthonormal parametrization. A third signed tie would have positive surface measure only if $`u_i,su_j,ru_k`$ were collinear. They are three distinct unit-sphere points, which cannot be collinear. All additional ties thus have zero interface measure. In particular the simultaneous zero projection at a planar base frame has codimension two and no interface mass.

For $`1<p<2`$, rotate coordinates smoothly in each bulk integral so that $`u_i\cdot x=x_1`$. Its possibly singular second derivative is bounded by a constant times
``` math
|x_1|^{p-2}(1+\|x\|^q)e^{-\|x\|^2/2}
```
for some finite $`q`$, which is integrable because $`p-2>-1`$. Derivatives of a path are bounded linear forms in these coordinates. Active indicators converge outside finitely many hyperplanes, so dominated convergence proves continuity of the bulk terms. One can first remove tubes around the zero hyperplanes in the differentiation; the integrals of the omitted singular terms tend to zero by this same bound. There is no zero-hyperplane delta mass, since the first derivative of $`|x_1|^p`$ is continuous at zero. For $`p\ge2`$ ordinary Gaussian-polynomial domination suffices.

The interface factor $`|z|^{p-1}`$ is continuous at zero for every $`p>1`$. Smooth interface parametrizations, bounded pair separations and Gaussian-polynomial domination give continuity of their integrals; active indicators converge off the additional null ties just discussed. Gaussian tails justify removing a preliminary large-radius truncation. The resulting bulk and surface Hessian expressions are continuous in local product-sphere coordinates, proving the lemma and the uniform Taylor remainder on bounded tangent sets. ◻

</div>

The argument deliberately does not cover $`p\le1`$: formal bulk differentiation then requires a different treatment of the zero hyperplanes.

# Calculation and strictness of the local quadratic form

## Angular, bulk and interface terms

The angular contribution in <a href="#eq:Q" data-reference-type="eqref" data-reference="eq:Q">[eq:Q]</a> follows by twice differentiating the gap expression, with gaps $`\alpha+t(a_{i+1}-a_i)`$. The first-order terms telescope. The second-order coefficient is $`H_{n,p}`$ as in <a href="#eq:HK" data-reference-type="eqref" data-reference="eq:HK">[eq:HK]</a>.

For a normal perturbation, the coefficient inside the winning cell of row $`i`$ is
``` math
\frac{p(p-1)}2|v_i\cdot x|^{p-2}(b_i\cdot x_\perp)^2
 -\frac p2\|b_i\|^2|v_i\cdot x|^p.
```
Integration over its two opposite sectors and the independent normal Gaussian gives the coefficient of $`\|b_i\|^2`$:
``` math
\begin{align*}
 D&=\frac{p m_{p-2}}\pi
 \left((p-1)\int_0^h\sin^{p-2}x\,dx
             -p\int_0^h\sin^p x\,dx\right)\\
 &=\frac{p m_{p-2}}\pi\sin^{p-1}h\cos h.
\end{align*}
```
We used $`m_p=p m_{p-2}`$ and integrated the derivative of $`\sin^{p-1}x\cos x`$, whose value at zero is zero for $`p>1`$.

Only neighbouring zero lines have active interfaces of positive measure at the cosine frame. On each such wall the common coordinate has magnitude $`\sin h\,|T|`$ and the row separation is $`2\cos h`$. The identity
``` math
\varphi_1(0)\mathbb E|T|^{p-1}=\frac{m_{p-2}}\pi
```
in Lemma <a href="#lem:C2" data-reference-type="ref" data-reference="lem:C2">3</a> gives coefficient $`-K_{n,p}`$. Consecutive rows use the sign $`s=-1`$; the wraparound uses $`s=+1`$. Hence the sum of interface terms is
``` math
-K_{n,p}\left(\sum_{i=0}^{n-2}\|b_i+b_{i+1}\|^2+
                  \|b_{n-1}-b_0\|^2\right)
 =-K_{n,p}\left(2\|B\|_F^2+\operatorname{tr}(B^{\mathsf T}A_nB)\right).
```
Since $`D=4K_{n,p}\cos^2h`$, combining bulk and surface terms proves the normal form. Mixed angular/normal terms vanish by oddness in the centered normal Gaussian coordinates, on both cells and interfaces. There is no linear normal term for the same reason. Lemma <a href="#lem:C2" data-reference-type="ref" data-reference="lem:C2">3</a> supplies the remainder.

## Spectrum and the orthogonal orbit

The antiperiodic Fourier vectors $`z_j=\exp(i(2k+1)\pi j/n)`$ diagonalize $`A_n`$, with eigenvalues $`2\cos((2k+1)\pi/n)`$. Its largest eigenvalue is $`2\cos\alpha`$, with multiplicity two and real eigenspace $`\operatorname{col}(V)`$. The next is $`2\cos3\alpha`$; for $`n=3`$ it is $`-2`$. The angular cycle Laplacian has kernel the constant vector and smallest positive eigenvalue $`4\sin^2\alpha`$. This proves positivity and the kernel assertion in Theorem <a href="#thm:main" data-reference-type="ref" data-reference="thm:main">1</a>. These zero velocities are precisely the common angular rotation and rotations of the base plane into its orthogonal complement.

For completeness, positivity transverse to that kernel yields actual local growth as follows. The orthogonal orbit is a compact embedded submanifold, even though its stabilizer $`O(d-2)`$ is nontrivial. Choose a closest orbit point to $`U`$ and rotate it to $`V`$. Put $`\Delta=U-V`$. Closest-point minimality implies that $`\Delta`$ is perpendicular to every orbit tangent. Its rowwise tangent projection $`w`$ onto the spheres remains perpendicular, since orbit tangents are rowwise tangent. Unit row lengths imply
``` math
\Delta_i\cdot v_i=-\tfrac12\|\Delta_i\|^2,\qquad
 \|w\|_F^2=D(U)+O(D(U)^2).
```
The chart used in <a href="#eq:path" data-reference-type="eqref" data-reference="eq:path">[eq:path]</a> has identity derivative on these tangent variables. Its Taylor expansion, stationarity and the transverse lower bound give
``` math
F_p(UU^{\mathsf T})-F_p(C_{\rm cos})
 \ge\lambda_{n,p}\|w\|_F^2+o(\|w\|_F^2).
```
Shrinking the neighbourhood gives <a href="#eq:growth" data-reference-type="eqref" data-reference="eq:growth">[eq:growth]</a>; invariance and compactness make it uniform along the orbit.

Finally take $`d=n`$. Positive semidefinite square roots $`C^{1/2}`$ depend continuously on $`C`$, also at a singular matrix. A factor of $`C_{\rm cos}`$ differs from the padded cosine frame by an orthogonal transformation. Thus all nearby correlation matrices have factors near the orbit, and the strict statement transfers to the full correlation domain. Discrete sign/permutation copies do not create extra infinitesimal kernels. This completes the proof of Theorem <a href="#thm:main" data-reference-type="ref" data-reference="thm:main">1</a>.

At $`p=2`$ the smaller transverse coefficient is the normal one,
``` math
\lambda_{n,2}=\frac{\tan(\alpha/2)}{\pi}(\cos\alpha-\cos3\alpha).
```
Its ratio to the angular coefficient is $`2\cos\alpha/(1+\cos\alpha)<1`$. Thus the gap weakens with increasing $`n`$; the result does not assert dimension-uniform stability. At $`n=4`$, $`A_4^2=2I`$ and $`A_4V=\sqrt2V`$, which gives <a href="#eq:four" data-reference-type="eqref" data-reference="eq:four">[eq:four]</a>.

# Four-coordinate comparisons

Write $`\mathcal F(C)=F_2(C)`$ and $`c_0=1-2\sqrt2/\pi`$ for four coordinates.

## A sharp independent-coordinate subclass

<div id="prop:block" class="proposition">

**Proposition 4**. *Suppose $`G_1,G_2,G_3`$ have a correlation matrix of rank at most two, and $`G_4`$ is an independent standard Gaussian. Then
``` math
\begin{equation}
\label{eq:block}
 \mathbb E\min_{1\leq i\leq4}G_i^2
 \geq 1-\frac{3}{2\pi}(4+\sqrt3-\sqrt{15})
 =0.11236067194906731320\ldots .
\end{equation}
```
Equality holds when the first three Gram lines are equally spaced. In particular, the bound is strictly larger than $`c_0`$.*

</div>

<div class="proof">

*Proof.* Conditioning on the independent coordinate preserves the planar tail comparison in Proposition <a href="#prop:planar" data-reference-type="ref" data-reference="prop:planar">2</a>. Thus three equally spaced directions minimise the objective in this entire class.

For that arrangement, write its squared minimum as $`R^2\sin^2\theta`$, with $`\theta`$ uniform on $`[0,\pi/6]`$, independent of $`R^2\sim\chi^2_2`$ and $`Z^2\sim\chi^2_1`$. For $`a\geq0`$,
``` math
\mathbb E\min(a^2R^2,Z^2)=2a^2\left(1-\frac{a}{\sqrt{1+a^2}}\right).
```
This follows by integrating the survival function of an exponential variable against that of $`Z^2`$. Therefore the value is
``` math
\frac{12}\pi\int_0^{\pi/6}
 \sin^2\theta\left(1-\frac{\sin\theta}{\sqrt{1+\sin^2\theta}}\right)d\theta.
```
Use $`u=\cos\theta`$ for the second term, and the antiderivative
``` math
\int\frac{1-u^2}{\sqrt{2-u^2}}du
 =\frac u2\sqrt{2-u^2}.
```
The result simplifies to <a href="#eq:block" data-reference-type="eqref" data-reference="eq:block">[eq:block]</a>. Its gap above $`c_0`$ is $`0.01267698810617338276\ldots`$. ◻

</div>

The same planar tail argument implies that a direct sum of two two-variable Gaussian blocks has objective at least that of four independent coordinates.

# Independence is a non-global local minimum

<div class="proposition">

**Proposition 5**. *For correlation matrices near $`I_4`$,
``` math
\begin{equation}
\label{eq:i-expansion}
 \mathcal F(C)=\mathcal F(I_4)+\kappa\sum_{i<j}c_{ij}^2+O(\|C-I_4\|^3),
 \qquad
 \kappa=\frac1\pi\int_0^\infty t e^{-t}
 \operatorname{erfc}(\sqrt{t/2})^2dt>0.
\end{equation}
```
Here $`\operatorname{erfc}(x)=(2/\sqrt\pi)\int_x^\infty e^{-s^2}ds`$. Thus $`I_4`$ is a strict local minimum, although
``` math
\mathcal F(I_4)=1-\frac{12-16/\sqrt3}{\pi}
 =0.12070214137740287227\ldots>c_0.
```*

</div>

<div class="proof">

*Proof.* Expand the joint Gaussian density around independence. The coefficient of $`c_{ij}`$ is odd in coordinates $`i,j`$ and integrates to zero against the minimum of the squares. Distinct-edge quadratic cross terms vanish by the same parity argument. The second derivative for a single edge is
``` math
\mathbb E\big[(Z_i^2-1)(Z_j^2-1)\min_k Z_k^2\big].
```
Use the layer-cake representation of the minimum and
``` math
\mathbb E[(Z^2-1)\mathbf1_{\{Z^2>t\}}]
 =2\sqrt t\,\varphi_1(\sqrt t).
```
The resulting second derivative is $`2\kappa`$. Analyticity of the density near $`I_4`$, with an integrable Gaussian domination for its derivatives, justifies the expansion. The value at $`I_4`$ also follows from the elementary evaluator below. ◻

</div>

This is an explicit obstruction to treating a successful local optimisation as a global certificate.

# An elementary evaluator for all four-variable frames

Assume first that distinct Gram lines are not parallel. For each pair $`(i,j)`$ and $`s=\pm1`$, set $`\rho=u_i\cdot u_j`$ and define
``` math
a=\sqrt{\frac{1+s\rho}{2}},\quad
 w=\frac{u_i+su_j}{\sqrt{2(1+s\rho)}},\quad
 n=\frac{u_i-su_j}{\sqrt{2(1-s\rho)}}.
```
Choose an orthonormal basis for $`\{n,w\}^{\perp}`$ and write, for the other two indices, $`\beta_k=u_k\cdot w`$ and $`b_k`$ for the coordinates of the residual projection. Let
``` math
W_{ij,s}=\frac{\mathbb E\big[|T|\,
 \mathbf1_{\{|\beta_kT+b_k\cdot Z|\geq a|T|\text{ for both }k\}}\big]}{\mathbb E|T|},
```
where $`T`$ and $`Z\in\mathbb R^{d-2}`$ are independent standard Gaussians. Then
``` math
\begin{equation}
\label{eq:interface-eval}
 \boxed{\mathcal F(UU^{\mathsf T})
 =1-\frac1\pi\sum_{i<j}\sqrt{1-\rho_{ij}^2}
 (W_{ij,+}+W_{ij,-}).}
\end{equation}
```
Indeed, the distributional Laplacian of the minimum has bulk part $`2`$ and an interface loss $`2|u_i\cdot x|\|u_i-su_j\|`$. Gaussian integration by parts and degree-two homogeneity give $`\mathbb E\Delta f=2\mathbb Ef`$. Integrating the interface loss gives <a href="#eq:interface-eval" data-reference-type="eqref" data-reference="eq:interface-eval">[eq:interface-eval]</a>.

Under the $`|T|`$ bias, the ratio $`r=Z/T`$ has density
``` math
\frac1{2(1+r^2)^{3/2}}\quad(d=3),\qquad
 \frac1{\pi(1+\|r\|^2)^2}\quad(d=4).
```
In dimension three, the forbidden conditions are two intervals. Their union probability is computed using the antiderivative $`r/(2\sqrt{1+r^2})`$.

In dimension four, each forbidden condition is a strip. Their intersection is a parallelogram $`K`$. The identity
``` math
\mathbb P(r\in K)=\frac1{2\pi}\int_{\partial K}\frac{r\cdot n}{1+\|r\|^2}\,ds
```
follows from the divergence theorem. For a counterclockwise edge starting at $`p`$, with unit tangent $`v`$ and length $`L`$, let $`h=\det(p,v)`$ and $`q=p\cdot v`$. Its contribution is
``` math
\frac{h}{2\pi\sqrt{1+h^2}}
 \left[\arctan\frac{q+L}{\sqrt{1+h^2}}
       -\arctan\frac{q}{\sqrt{1+h^2}}\right].
```
Inclusion–exclusion for the two strips completes a finite elementary evaluation. Singular configurations follow by taking limits; duplicated lines can be removed. The implementation uses floating-point rank thresholds, not symbolic limit evaluation. Its numerical output is not an outward-rounded certificate.

# A geometric reduction and independent evaluator

Let $`K_S=\operatorname{conv}\{\pm u_i:i\in S\}`$. Write $`V_2(K)`$ for the second intrinsic volume, normalised to be ordinary area in dimension two and half the surface area in dimension three. Long’s known identity specializes, for a polytope whose vertices have norm one, to
``` math
\begin{equation}
\label{eq:intrinsic}
 \mathbb E\max_i(u_i\cdot Z)^2=1+\frac{V_2(K)}\pi.
\end{equation}
```
For completeness, its support function satisfies $`h^2+\|\nabla_Sh\|^2=1`$ almost everywhere on the sphere. The second intrinsic-volume formula is
``` math
V_2(K)=\frac1{2\kappa_{d-2}}
 \int_{S^{d-1}}\big((d-1)h^2-\|\nabla_Sh\|^2\big)\,d\sigma.
```
Here $`\kappa_m`$ is the volume of the unit $`m`$-ball and $`d\sigma`$ is unnormalised surface measure. Substitution, $`|S^{d-1}|=2\pi\kappa_{d-2}`$, and Gaussian radial integration prove <a href="#eq:intrinsic" data-reference-type="eqref" data-reference="eq:intrinsic">[eq:intrinsic]</a>. The intrinsic-volume formula can be obtained by expanding the parallel-body volume to its quadratic support-function coefficient and integrating by parts on the sphere; approximation extends it from smooth support functions to polytopes.

The pointwise inclusion–exclusion identity for a minimum gives
``` math
\begin{equation}
\label{eq:ie-v2}
 \mathcal F(C)=1+\frac1\pi\left(-\sum_{|S|=2}V_2(K_S)
                  +\sum_{|S|=3}V_2(K_S)-V_2(K_{[4]})\right).
\end{equation}
```
For a pair, $`V_2(K_{ij})=2\sqrt{1-c_{ij}^2}`$. For a triple $`T`$, choose four sign patterns with its first sign positive. Define
``` math
R_{T,\varepsilon}=
 \big\|(\varepsilon_j u_j-\varepsilon_i u_i)
       \wedge(\varepsilon_k u_k-\varepsilon_i u_i)\big\|,
```
where the wedge norm is twice the triangle area. Then $`V_2(K_T)=\frac12\sum_{\varepsilon}R_{T,\varepsilon}`$, also by continuity in degenerate triples.

## The residual rank-three inequality

In rank three, $`\mathcal F(C)\ge c_0`$ is equivalent to
``` math
\begin{equation}
\label{eq:rank3-residual}
 4\sum_{i<j}\sqrt{1-c_{ij}^2}
 -\sum_{T,\varepsilon}R_{T,\varepsilon}
 +\operatorname{area}(\partial K_{[4]})\leq4\sqrt2.
\end{equation}
```
For a generic configuration, $`K_{[4]}`$ has twelve triangular faces, paired antipodally. Of the sixteen unsigned triangle patterns, six are its face pairs and ten are not. Hence <a href="#eq:rank3-residual" data-reference-type="eqref" data-reference="eq:rank3-residual">[eq:rank3-residual]</a> becomes a radical inequality involving the ten non-face patterns:
``` math
\sum_{\text{ten non-face patterns}}R_{T,\varepsilon}
 \geq4\sum_{i<j}\sqrt{1-c_{ij}^2}-4\sqrt2.
```
For a triple with signed vertices $`v_i,v_j,v_k`$, set $`N=(v_j-v_i)\times(v_k-v_i)`$. With $`l`$ the omitted index, it is a non-face pattern precisely when $`|u_l\cdot N|>|v_i\cdot N|`$. At four-way ties, the generic ten-pattern rule must be replaced by the surface-area formula or a continuous limit. In particular, simply including every equality case in this indicator is not valid.

## The residual rank-four inequality

If $`U`$ is invertible, extend each triple sign pattern by setting the omitted sign to $`+1`$ and $`-1`$, giving $`\varepsilon^+`$ and $`\varepsilon^-`$. The adjacent facet normals are proportional to $`U^{-1}\varepsilon^+`$ and $`U^{-1}\varepsilon^-`$. Let their angle be $`\alpha_{T,\varepsilon}`$. The normal cone has two-dimensional Gaussian angle $`\alpha/(2\pi)`$, so
``` math
V_2(K_{[4]})=\sum_{T,\varepsilon}
 R_{T,\varepsilon}\frac{\alpha_{T,\varepsilon}}{2\pi}.
```
Consequently $`\mathcal F(C)\ge c_0`$ is equivalent to
``` math
\begin{equation}
\label{eq:rank4-residual}
 \sum_{T,\varepsilon}R_{T,\varepsilon}
       \left(1-\frac{\alpha_{T,\varepsilon}}\pi\right)
 \geq4\sum_{i<j}\sqrt{1-c_{ij}^2}-4\sqrt2.
\end{equation}
```
Equations <a href="#eq:rank3-residual" data-reference-type="eqref" data-reference="eq:rank3-residual">[eq:rank3-residual]</a> and <a href="#eq:rank4-residual" data-reference-type="eqref" data-reference="eq:rank4-residual">[eq:rank4-residual]</a> are explicit finite-dimensional proof obligations. **Neither is proved globally in this note.** The geometric implementation uses a convex-hull surface-area computation in rank three, so it also handles non-triangular faces.

# An independent-noise exclusion region

<div class="lemma">

**Lemma 6** (Independent-noise floor). *If $`C-D`$ is positive semidefinite for $`D=\operatorname{diag}(d_1,\ldots,d_n)`$, $`d_i\ge0`$, then for $`G\sim N(0,C)`$ and independent standard normals $`Z_i`$,
``` math
\mathbb E\min_i|G_i|^p\ge\mathbb E\min_i|\sqrt{d_i}Z_i|^p,\qquad p>0.
```*

</div>

<div class="proof">

*Proof.* Write $`G=Y+D^{1/2}Z`$ with independent Gaussian vectors. An even decreasing Gaussian density places the most mass in an interval of fixed length when the interval is centered. Hence $`\mathbb P(|\sqrt dZ+m|>a)\ge\mathbb P(|\sqrt dZ|>a)`$ for $`a\ge0`$, any shift $`m`$, and $`d>0`$. Conditional on $`Y`$, multiply these inequalities using independence. Average and integrate against $`pa^{p-1}da`$. If any $`d_i=0`$, the right side is zero and the conclusion immediate. ◻

</div>

This elementary supporting comparison is not claimed as a new general principle. For four unit variances it gives $`\mathcal F(C)\ge\lambda_{\min}(C)\mathcal F(I_4)`$. Thus the target holds whenever
``` math
\lambda_{\min}(C)\ge
 \tau:=\frac{1-2\sqrt2/\pi}{1-(12-16/\sqrt3)/\pi}.
```
The sufficient condition $`\max_i\sum_{j\ne i}|c_{ij}|\le1-\tau`$ follows from the Gershgorin bound. This excludes an interior region, not the entire remaining domain.

# Computational evidence and limits

The universal proofs do not depend on numerical optimization. The evaluators implement four-row frames, not the all-$`n`$ theorem. The original experiment compared them on 100 sampled frames of each rank three and four, with maximum differences $`9.77\cdot10^{-15}`$ and $`6.00\cdot10^{-15}`$. These are floating-point consistency tests, not error enclosures. The original interface evaluator could return non-finite results on valid near-singular inputs. The repaired version raises `NumericalDomainError` on unsupported duplicate or near-singular frames; it does not claim full numerical coverage. Finite-precision rank thresholds are not symbolic limiting arguments. Explicit finiteness checks precede tolerance comparisons, and missing or malformed saved batches are rejected.

The original BFGS experiment used 100 starts in each dimension, sampled from an angular box rather than uniformly over correlations. Endpoints were classified by objective value within $`10^{-6}`$, not by geometric basin or certified stationarity:

<div class="center">

| Dimension | Cosine | Triple-plus-one | Independence | Success flags |
|:----------|-------:|----------------:|-------------:|--------------:|
| Three     |     28 |              72 |            0 |       100/100 |
| Four      |     31 |              58 |           11 |        78/100 |

</div>

All 22 unsuccessful exits are retained. The two lowest values were $`0.09968368384294646`$ and $`0.09968368384290127`$. Reevaluation of saved endpoints is distinct from regenerating optimization trajectories. These samples cannot supply a probability of global correctness.

Eight original finite-difference paths for <a href="#eq:four" data-reference-type="eqref" data-reference="eq:four">[eq:four]</a> gave a largest coefficient discrepancy $`3.16\cdot10^{-6}`$ at step $`0.002`$. This does not certify the uniform remainder. Exact algebra over $`\mathbb Q(\sqrt2)`$ checks the four-cycle identities and rejects an altered matrix; it does not check analytic differentiation. The source release supplies computational inputs, commands and execution records separately from the analytic proof.

The repaired scripts were replayed in normal and optimized Python modes, including all 200 saved endpoints and seven guard-test methods. Checks use explicit exceptions or test-framework checks rather than optimization-removable assertions. The runtime was Python 3.11.14, NumPy 2.4.6, SciPy 1.17.1 and SymPy 1.14.0, differing from the original numerical pins; this is cross-version validation, not exact historical-environment reproduction. The original optimization trajectories were not rerun. Historical termination messages and final gradients were not recorded and are not reconstructed. Deterministic difficult-frame tests document explicit refusal by the elementary evaluator and finite numerical comparison by the intrinsic evaluator, not rigorous error enclosures.

# Remaining questions

The local theorem has an existential radius, not an explicit certified radius. No exhaustive interval cover of the correlation domain is supplied. Independence is itself a nonglobal strict local minimum, so every nontrivial straight covariance ray from it initially moves uphill. This prevents a naive monotone rank-reduction argument, although it does not exclude other endpoint comparisons.

The rank-three nonface choices and rank-four angle weights remain coupled geometric quantities. Their exact formulas are obligations, not proofs of the required inequalities. The noise-floor lemma excludes a specified interior region only. The global four-coordinate bound and Kunisky’s stochastic conjecture remain unresolved here. Orders $`0<p\le1`$ also require a different treatment of zero-hyperplane singularities. The literature comparisons do not justify first-ever priority claims.

# Availability and declarations

This manuscript and the accompanying source release are available under CC0. The mathematical arguments do not depend on access to private review records. AI assistance was used for mathematical development, checking and exposition; this is not proof-assistant certification or external journal peer review. Attribution is Anonymous. The work involves no human participants, personal data or animal experiments. Author-specific funding, contribution and competing-interest information is not supplied; none is inferred.

<div class="thebibliography">

9 D. Kunisky. *A revision of Litvak’s conjecture on Gaussian minima and a volumetric zone conjecture*. arXiv:2605.02023v1, 2026. <https://arxiv.org/abs/2605.02023v1>. Y. Long. Hyperbolic embedding of infinite-dimensional convex bodies. *Geometriae Dedicata* 219, article 17, 2025. <https://doi.org/10.1007/s10711-024-00957-5>. Inspected version: <https://arxiv.org/html/2305.13428v3>. A. E. Litvak. Around the simplex mean width conjecture. In *Analytic Aspects of Convexity*, Springer INdAM Series 25, 73–84, 2018. <https://doi.org/10.1007/978-3-319-71834-7_5>. Y.-F. Zhong, T. Ma and Z.-C. Hu. *A note on the Gaussian minimum conjecture*. arXiv:2008.06211v1, 2020. <https://arxiv.org/abs/2008.06211v1>. cyclemath (editor and releaser). *Convex Hulls of Great Subspheres and Gaussian Minima*. Preliminary AI-assisted manuscript, public version 0.1, September 2026. <https://github.com/cyclemath/great-subspheres-gaussian-minima>.

</div>
