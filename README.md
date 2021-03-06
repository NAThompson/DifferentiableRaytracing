slidenumbers: true

## Raytracing Differentiable Surfaces

![inline](figures/ornl_logo.svg)

Nick Thompson

---

![](figures/albrecht_durer.jpg)

---

## How to render a sphere

First: Compute the intersection of a ray and a sphere.

Ray is $$\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$$, sphere: $$(\mathbf{P} - \mathbf{c})^2 = R^2$$.

$$
\mathbf{d}^2t^2 + 2(\mathbf{o} - \mathbf{c})\cdot \mathbf{d}t + (\mathbf{o}-\mathbf{c})^2 = R^2
$$

A quadratic for $$t$$. No real roots $$\implies$$ ray doesn't hit sphere.

---

![original](figures/shirley_book1.jpg)

[.footer: Final scene from Peter Shirley's *Raytracing in a Weekend* series. You can do a lot with a good random number generator, Snell's law, and the ability to solve quadratics.]

---

## Rendering Ellipsoids

$$\frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} = 1$$

A ray $$\mathbf{r}(t) = \mathbf{o} + t\mathrm{d}$$ intersects this at

$$
\frac{(o_x+td_x)^2}{a^2} + \frac{(o_y+td_y)^2}{b^2} + \frac{(o_z+td_z)^2}{c^2} = 1
$$

Again, just a quadratic.

---

![inline](figures/ellipsoid.png)

---

![original](figures/ellipsoid_curvature.gif)


[.footer: *An ellipsoid pseudocolored by Gaussian curvature*]

---

## Implicit Equation of a Torus

$$(x^2 + y^2 + z^2 + R^2 - r^2)^2 = 4R^2(x^2+y^2)$$

Plug in $$x = o_x + td_x, y = o_y + td_y, z = o_z + td_z$$ and expand to get a quartic equation for $$t$$.

- No real roots implies no intersection of ray with torus
- Multiple real roots implies multiple intersections; take the smallest positive real root.

---

![left](figures/torus_curvature.png)

Gaussian curvature of torus computed from numerically evaluated first and second fundamental forms. Rendered via solution of quartic equation.

---

![left](figures/pv_torus_curvature.png)

Gaussian curvature of torus computed using the curvature filter in Paraview.
Artifacts from the triangle mesh are clearly visible and the curvature is incorrect everywhere.

---

## Raytracing *differentiable* objects

Naive raytracing is *agnostic* to the differentiability of the objects in your scene.

All that matters to raytracing is that you can compute the intersection of your objects with a line.

---

## Ok, why do we need to focus on raytracing *differentiable* structures?

Raytracing triangle meshes has become dominant in game programming. The assumption is that *creation* of the mesh is not that expensive-a good assumption for game asset generation; a bad one for visualizing ODE/PDE solutions.

CAD programs render diagrams using smoother surfaces; e.g., Bezier and Catmull-Clark.

---

## A billion triangles a frame!

![inline](https://youtu.be/qC5KtatMcUw?t=106)

---

## Benefits of exploiting smoothness

- Memory requirements are drastically reduced relative to triangle meshes.
- No error prone conversion into triangle meshes.
- Pseudocoloring by derived quantities results in fewer visual artifacts.
- No need to subsample high-accuracy numerical solutions

^ High accuracy numerical methods might produce sparse pointsets in space. Without smooth interpolation, we have to subsample to extreme density to avoid artifacts.

---


## Implicit equations are "easy" to render

Well, relatively easy, since 1D rootfinding problems are easy. But many interesting objects, such as Beziers and B-splines are parametric surfaces. Rendering them is a multivariate rootfinding problem-much harder.

Let's focus on those.


---

Parametric equation of helicoid

$$\sigma(u, v) =
\begin{pmatrix}
rv\cos(2\pi u)\\
rv\sin(2\pi u)\\
\lambda (u-1/2)
\end{pmatrix}$$

Intersects ray when

$$\mathbf{o} + t\mathbf{d} = \sigma(u,v)$$

^ We must recover $$(u,v,t)$$ to render, in other words we must solve a multivariate rootfinding problem.

![left](figures/helicoid_smooth_cool_warm_lambertian.png)


---

## Helicoid ray intersections

![left](figures/helicoid_bounded.png)

Any method we use to compute the parameters $$(u,v,t)$$ will require a good initial guess and bounds on the parameters.

Let's surround the helicoid with a cylinder $$x^2 + y^2 = r^2$$.

Intersection points on the cylinder will give us bounds on $$u,v$$ and $$t$$.

---

## Inequalities are just as important as numerics!

Working to *constrain* a ray intersection with an object is just as important as the convergence of the method to the correct solution.

The ability to say "this ray *does not* intersect our object" in as few flops as possible is just as important for your performance as the numerical method to compute the intersection when it does.

---

What if we *can't* constrain our initial guess very tightly?

The Newton fractal teaches us that predicting *which* root Newton's method will converge to is difficult.

In ray tracing, we need the *minimal* root to determine the closest intersection point to the ray.

![left](figures/newton_fractal_z8_15z4_m16_viridis.png)

[.footer: *Newton fractal for $$z^8 - 15z^4 - 16$$*]

---



![left](figures/wrong_solutions.png)

A helicoid pseudocolored by Gaussian curvature. The pseudocoloring changes discontinuously where the non-minial roots are identified.

---

## Multivariate root-finding is the worst!

![inline](figures/the_worst.gif)


---

## Multivariate root-finding is the worst!

Bounced rays make the problem even worse than the generic multivariate rootfinding problem.

A ray that bounces off an object satisfies the equation $$f(\mathbf{o} + t\mathbf{d}) = 0$$ at $$t = t_{\min}$$, and it *might* hit the object again.

So you need to find the $$t$$ *closest* to $$t_{\min}$$ that hits the object, but not the $$t = t_{\min}$$ solution.

---

And numerical imprecision might mean your $$t = t_{\min} + \epsilon$$ solution is essentially the *same* as your $$t = t_{\min}$$ solution.

This generates a visual artifact known as [shadow acne](https://digitalrune.github.io/DigitalRune-Documentation/html/3f4d959e-9c98-4a97-8d85-7a73c26145d7.htm).

![left](figures/shadow_acne.jpg)


[.footer: *Image courtesy of Digital Rune.*]

---

## Let's examine three methods to render parametric surfaces

- Newton's method with backtracking.
- The multivariate Halley iterate.
- A second order geometric method of Ko et. al.


---

## Newton's method

Given $$\mathbf{f}(u,v,t) := \sigma(u,v) - \mathbf{o} - t\mathbf{d}$$, the Newton's method is

$$(u_{k+1}, v_{k+1}, t_{k+1}) = (u_{k}, v_{k}, t_{k}) -\mathbf{J}_{\mathbf{f}}^{-1}\mathbf{f}(u_k,v_k,t_k)$$

---

## Newton's method is bad

- Need to compute derivatives/Jacobians
- Slow with multiple roots
- Full-measure regions of initial conditions where it diverges.
- Can converge to values which are not roots.
- For ray intersections, it *ignores* your initial guess for $$t$$.

^ Since Newton's method computes the intersection of line with tangent plane, $$t$$ is uniquely determined by the guesses for the parametric coorindates $$(u_0, v_0)$$.

But let's examine a hack that makes it a bit better.

---

## "Patched" Newton's method: NR§9.7

Given $$\mathbf{f}\colon \Omega \subset \mathbb{R}^3 \to \mathbb{R}^3$$ defined as $$\mathbf{f}(u,v,t) = \sigma(u,v) - (\mathbf{o} + t\mathbf{d})$$, we wish to find $$(u,v, t)$$ s.t. $$\mathbf{f}(u,v,t) = \mathbf{0}$$.

Newton update: $$(u_{k+1}, v_{k+1}, t_{k+1}) = (u_{k}, v_{k}, t_{k}) -\mathbf{J}_{\mathbf{f}}^{-1}\mathbf{f}(u_k,v_k,t_k)$$.

Patch: Define $$f(u,v,t) := \frac{1}{2}\left<\mathbf{f}(u,v,t), \mathbf{f}(u,v,t) \right>$$.

---

## Patched Newton's method

If $$f(u_{k+1}, v_{k+1}, t_{k+1}) > f(u_{k}, v_{k}, t_k)$$ or $$(u_{k+1}, v_{k+1}, t_{k+1}) \not \in \Omega$$, reject the Newton step.

Define $$g(\lambda) := f(u_{k} + \lambda \delta u, v_{k} + \lambda \delta v, t_{k} + \lambda \delta t)$$, $$\lambda \in (0,1)$$.

Choose $$\lambda$$ that minimizes $$g$$.

---

![left](figures/helicoid_newton_wbacktracking.png)

Helicoid rendered with Newton's method patched with backtracking. Initial guess for $$t$$ provided by bounding cylinder intersection, no bounced rays.

^ Interval arithmetic can be used to determine existence and convergence of Newton's method. See Toth.
^ This is based on ideas from Krawczyk's operator, which bounds volumes of convergent Newton iterates.
^ If subsequent volumes do not intersect, there is no solution!

---

![left](figures/helicoid_newton_wbacktracking.png)
![right](figures/helicoid_newton_cost_inferno.png)



[.footer: Helicoid rendered via Newton's method juxtapositioned against pixel cost map. Lighter colors are more expensive.]

^ When we do not intersect the cylinder, we get very cheap pixels. When we hit the helicoid, but no solution exists, we get very expensive pixels.

---

## Performance

Since we are restricted to $$n \le 4$$ in computer graphics, we can soberly consider rootfinding methods that are asymptotically ill-advised.

e.g., multivariate Halley's method, which requires formation of Hessian tensors.

---

## Halley's method

$$t_{k+1} = t_{k} - \frac{2f(t_k)f'(t_k)}{2f'(t_k)^2 - f(t_k)f''(t_k)}$$

Convergence is cubic, much more importantly, *which* root Halley's method converges to is much more predictable; see [here](http://faculty.nps.edu/bneta/papers/BasinsThirdFourth.pdf) and [here](http://dx.doi.org/10.1016/j.amc.2011.07.076).

---

![left](figures/newton_fractal_viridis_z5m1.png)
![right](figures/halley_fractal_viridis_z5m1.png)


[.footer: *Newton (left) vs Halley fractal for $$z^5 -1$$.*]

---

## Multivariate Halley Iterate

Given $$\mathbf{f}\colon \Omega \subset \mathbb{R}^{3} \to \mathbb{R}^3$$, (e.g., $$\mathbf{f}(u,v,t) = \sigma(u,v) - \mathbf{o} - t\mathbf{d}$$), define

Then the Jacobian $$\mathbf{J}_f \colon \Omega \subset \mathbb{R}^{3} \to \mathbb{R}^{3\times 3}$$ is a linear operator for each $$\mathbf{q} \in \Omega$$.

And the Hessian $$\mathbf{H}_f\colon \Omega \subset \mathbb{R}^{3} \to \mathbb{R}^{3\times 3 \times 3}$$ is a bilinear operator for each $$\mathbf{q} \in \Omega$$.

---

## Multivariate Halley iterate

Define $$\mathbf{a}^{\nu} \in \mathbb{R}^3$$ as the solution to $$\mathbf{J}_{f(\mathbf{q}^{\nu})}\mathbf{a}^{\nu} = -\mathbf{f}(\mathbf{q}^{\nu})$$.

Define $$\mathbf{b}^{\nu} \in \mathbb{R}^3$$ as the solution to  $$\mathbf{J}_{f(\mathbf{q}^{\nu})}\mathbf{b}^{\nu} = \mathbf{H}_{f(\mathbf{q}^{\nu})}(\mathbf{a}^{\nu}, \mathbf{a}^{\nu})$$.

Then the multivariate Halley iterate is

$$\mathbf{q}^{\nu + 1} = \mathbf{q}^{\nu} + (\mathbf{a}^{\nu} \otimes \mathbf{a}^{\nu})\oslash(\mathbf{a}^{\nu} + \mathbf{b}^{\nu}/2)$$.

The $$\otimes$$ and $$\oslash$$ are defined componentwise; see [Cuyt](https://dl.acm.org/doi/pdf/10.1145/3147.3162).

---

![left](figures/halley_helicoid.png)

Helicoid rendered via the multivariate Halley iterate.

Speed is 1.5x of the Newton iterate.

---

![left](figures/halley_helicoid.png)
![right](figures/helicoid_halley_cost_inferno.png)


---

![left](figures/helicoid_newton_cost_inferno.png)
![right](figures/helicoid_halley_cost_inferno.png)

[.footer: *Pixel cost map for Newton's method (left) vs Halley's method.*]

---

## Is it the *method*, or the *problem*?

Both the Halley and Newton method struggle when the Jacobian is nearly singular.

Is this a fundamental aspect of surface intersection, or incidental?

---

## Conditioning of multivariate rootfinding

Let $$w^{*}$$ be a root of $$\mathbf{f}$$, and let $$w$$ be a root of a function $$\tilde{\mathbf{f}}  = \mathbf{f} + \Delta \mathbf{f}$$. Then

$$\frac{\left\| w - w^{*} \right\|}{ \left\| w^{*} \right\| } \le \frac{\left\| J_{\mathbf{f}}(w^{*})^{-1} \right\| }{ \left\| w^{*} \right\| }\left\|\mathbf{f}(w) \right\| =: \mathrm{cond}(\mathbf{f},w^{*}) \left\|\mathbf{f}(w) \right\|$$

See [Corless](https://link.springer.com/book/10.1007/978-1-4614-8453-0). (Aside: this also solves the problem of shadow acne.)


---

## In English

If you're surveying the Bonneville salt flats, you don't blame Newton's method when your laser hits the Wasatch front.

![](figures/bonneville_salt.jpg)

---

![left](figures/helicoid_condition_number.png)

A helicoid pseudocolored by the logarithm of the condition number of the multivariate rootfinding problem

$$\mathbf{f}(u,v,t) := \sigma(u,v) - \mathbf{o} - t\mathbf{d} = 0$$

---

![left](figures/helicoid_condition_number.png)
![right](figures/helicoid_halley_cost_inferno.png)

[.footer: Helicoid pseudocolored by logarithm of condition number juxtapositioned against the pixel cost map for the Halley iterate.]
---

The inverse of the Jacobian at the root constrains our accuracy *whether we use it explicitly or not*.

A *better criticism* of the Newton or Halley iterate is that they can fail because a singular Jacobian is encountered at any point *away* from the root.

---

## Time for New Ideas

Both the Newton and Halley iterate are *generic* multivariate rootfinding algorithms.

Can we find a new method that explicitly exploits the geometry of the ray intersection problem?

---

![inline](figures/buckle_up.gif)

## It's time for your daily dose of differential geometry

---

## A geometric surface intersection algorithm from [Ko et. al.](https://doi.org/10.1016/j.cagd.2013.07.001)

![left](figures/ko_method_graphic_problem.png)

Where does the orange ray intersect the surface?

---

![left](figures/ko_method_initial_guess.png)

Start with an initial guess $$\mathbf{p}_0 = \sigma(u_0, v_0)$$.

---

![left](figures/ko_method_ray_plane.png)

Compute the plane containing the ray $$\mathbf{o} + t \mathbf{d}$$ and $$\mathbf{p}_0 = \sigma(u_0, v_0)$$.

Its normal is $$(\mathbf{p}_0-\mathbf{o})\times \mathbf{d}$$.

The solution lies on this plane, but we don't know where!

We want to construct a curve $$\gamma$$ that lies in the plane and the surface.

---

## Formalizing

Given a parametric surface $$\sigma \colon \Omega \subset \mathbb{R}^2 \to \mathbb{R}^3$$ and an initial guess $$(u_0, v_0)$$ for an intersection with a ray $$\mathbf{o} + t\mathbf{d}$$, consider a unit speed curve

$$\gamma\colon (-\delta, \delta) \to \sigma(\Omega) \cap \mathcal{P}, \gamma(0) = \sigma(u_0, v_0)$$

where $$\mathcal{P}$$ is the plane containing the ray and $$\sigma(u_0,v_0)$$.

---

## Taylor expand

$$
\gamma(s) = \gamma(0) + s\dot{\gamma}(0) + \frac{s^2}{2}\ddot{\gamma}(0) + \mathcal{O}(s^3) = \sigma(u_0, v_0) + s\hat{\mathbf{t}} + \frac{\kappa}{2}s^2 \hat{\mathbf{n}} + \mathcal{O}(s^3)
$$

How does the surface and the plane influence $$\dot{\gamma}(0)$$ and $$\ddot{\gamma}(0)$$?

---

Let $$\mathbf{n}_{\mathcal{P}}$$ be the normal to the plane containing the ray and $$\sigma(u_0, v_0)$$, and let $$\mathbf{n}_{\sigma}$$ be the normal to $$\sigma$$ at $$(u_0,v_0)$$.

$$\mathbf{n}_{\mathcal{P}} \propto (\mathbf{o} - \sigma(u_0,v_0)) \times \mathbf{d}$$.

$$\mathbf{n}_{\sigma} \propto \partial_u \sigma \times \partial_v \sigma$$.

$$\dot{\gamma}(0)$$ must be orthogonal to both of them, so take $$\dot{\gamma}(0) \propto \mathbf{n}_{\mathcal{P}}\times \mathbf{n}_{\sigma}$$, and normalize.

---

$$\hat{\mathbf{t}} := \dot{\gamma}(0) \in \mathcal{T}_{\gamma(0)}\sigma$$, so there exist constants $$\dot{u}, \dot{v} \in \mathbb{R}$$ such that $$\hat{\mathbf{t}} = \dot{u}\partial_{u}\sigma + \dot{v}\partial_{v}\sigma$$. Dot with $$\partial_{u}\sigma$$, then $$\partial_v\sigma$$ to obtain

$$\begin{pmatrix}
\hat{\mathbf{t}}\cdot \partial_u\sigma \\
\hat{\mathbf{t}}\cdot \partial_v\sigma
\end{pmatrix}=\begin{pmatrix}
E & F \\
F & G
\end{pmatrix}\begin{pmatrix}
\dot{u} \\
\dot{v}
\end{pmatrix}$$

where $$E := \left\|\partial_u\sigma\right\|^2, F := \left<\partial_u\sigma, \partial_v\sigma \right>, G:= \left\|\partial_v\sigma\right\|^2$$ are the first fundamental forms.

*This matrix only fails to be invertible when the surface is irregular*.

---

## What do we know about $$\ddot{\gamma}(0)$$?

Since $$\dot{\gamma}$$ is unit speed, $$\frac{d}{dt}\dot{\gamma}(t)^2 = 2\dot{\gamma}(t)\cdot \ddot{\gamma}(t) = 0$$, so it is orthogonal to the tangent.

Moreover, $$\gamma$$ lies in a plane, so its torsion $$\tau$$ must vanish. Its binormal therefore must be $$\pm \mathbf{n}_{\mathcal{P}}$$.

So $$\ddot{\gamma}(0) \perp \dot{\gamma}(0)$$ and $$\ddot{\gamma}(0) \perp \mathbf{n}_{\mathcal{P}}$$.

---

So $$\ddot{\gamma}(0) \propto \mathbf{n}_{\mathcal{P}} \times \dot{\gamma}(0)$$. But what is the constant of proportionality?

Let's just pull this one out of a differential geometry book . . .


---

$$\kappa\hat{\mathbf{n}} = \frac{\kappa_{\sigma}}{1- (\mathbf{n}_{\mathcal{P}}\cdot \mathbf{n}_{\sigma})^2}(\mathbf{n}_{\sigma} -(\mathbf{n}_{\mathcal{P}}\cdot \mathbf{n}_{\sigma}) \mathbf{n}_{\mathcal{P}})$$


$$\kappa_{\sigma} := e\dot{u}^2 + 2f\dot{u}\dot{v} + g\dot{v}^2$$,

where

$$e:= \hat{\mathbf{n}}_{\sigma} \cdot\partial_{u}^2\sigma, f := \hat{\mathbf{n}}_{\sigma} \cdot\partial_{u}\partial_{v}\sigma, g:= \hat{\mathbf{n}}_{\sigma} \cdot\partial_{v}^2\sigma$$

are the [second fundamental forms](https://en.wikipedia.org/wiki/Second_fundamental_form).

---

### Iteration: Intersect curve with line

$$\gamma(0) + s\hat{\mathbf{t}} + \frac{\kappa}{2}s^2 \hat{\mathbf{n}} = \mathbf{o} + t\mathbf{d}$$

Dot with $$\hat{\mathbf{t}}$$ and $$\mathbf{d}$$ and eliminate $$t$$ to obtain

$$
\frac{\kappa}{2}s^2 - \frac{\mathbf{d}\cdot\hat{\mathbf{n}}}{\mathbf{d}\cdot\hat{\mathbf{t}}} s = (\mathbf{o} - \gamma(0))\cdot \hat{\mathbf{n}} -  \frac{\mathbf{d}\cdot\hat{\mathbf{n}}}{\mathbf{d}\cdot\hat{\mathbf{t}}}(\mathbf{o} -\gamma(0))\cdot\hat{\mathbf{t}}
$$

Choose root than minimizes $$t$$.

---

### Iteration: Intersection curve with line

The update is

$$u_{k+1} := u_{k} + s\dot{u}, v_{k+1} := v_{k} + s\dot{v}$$.

---

![left](figures/helicoid_ko_method.png)

A helicoid rendered via the geometric method of Ko et. al.

No Jacobians are inverted in this calculation, yet the non-minimal $$t$$ solutions are found where the Jacobian is close to singular.

---

![left](figures/helicoid_ko_method.png)
![right](figures/helicoid_condition_number.png)

[.footer: A well-conditioned root *behind* a poorly-conditioned root causes the worst problems.]

---

## Aside: Residuals

Raytracers are *brutal* testers of numerical codes, and not in a way you always want.

Test residuals of every intersection!

---

## Implicit surface intersection residuals

For an intersection $$(\hat{x}, \hat{y}, \hat{z})$$ on an implicit surface $$f(x,y,z) = 0$$, we have an expected residual
$$f(\hat{x}, \hat{y}, \hat{z}) = f(x(1+\epsilon_x), y(1+\epsilon_y), z(1+\epsilon_z)) = x\epsilon_x \partial_x f + y\epsilon_y \partial_y f + z\epsilon_z \partial_z f + \mathcal{O}(\left\|\epsilon\right\|^2)$$

So a reasonable check is

$$
\left| f(\hat{x}, \hat{y}, \hat{z}) \right| \le \mu\left(|\hat{x}\partial_x f| + |\hat{y}\partial_y f| + |\hat{z}\partial_z f| \right)
$$

---

## Parametric surface intersection residuals

For the intersection of a ray $$\mathbf{o} + t\mathbf{d}$$ with a parametric surface $$\sigma(u,v)$$, we recover $$(\hat{t}, \hat{u}, \hat{v})$$.

The intersection point can be computed two ways: $$\hat{\mathbf{p}}_1 = \mathbf{o} + \hat{t}\mathbf{d}$$ or $$\hat{\mathbf{p}}_2 = \sigma(\hat{u}, \hat{v})$$.

Compute it *both* ways to get a residual, and verify that

$$
\left\|\hat{p}_1 - \hat{p}_2  \right\| \le \mu \left(|\hat{t}| \left\|\mathbf{d}\right\| + |\hat{u}| \left\| \partial_u \sigma \right\| + |\hat{v}| \left\| \partial_v \sigma \right\| \right)
$$

---

## Vision

Without clean abstractions, we cannot form teams of teams.

These techniques are a technical path to ensuring users of SciViz have only a single abstraction they are responsible for: Constructing a method of evaluating a continuous function.

---

## Vision: The *only* technical knowledge other teams should need to have to use our software

I have $$f\colon \mathbb{R}^3 \to \mathbb{R}$$. I can isocontour and volume render. Pseudocode:

```cpp
auto f = [&](double x, double y, double z) -> double {
    // implementation or binding to your function
    };
auto actor1 = vtkm::isocontour(f, domain, iso_value, pseudocolor_lambda);
auto actor2 = vtkm::volume_render(f, domain, another_pseudocoloring);
```

![](figures/ProtonVis.png)

---

## Vision

I have a parametric surface $$\sigma \colon \mathbb{R}^{2} \to \mathbb{R}^3$$. I can render geometry. Pseudocode:

```cpp
auto pseudocolor = [&](double u, double v) { return viridis(sigma.gaussian_curvature(u,v)); };
auto actor = vtkm::parametric_surface(sigma, pseudocolor);
auto image = vtkm::render(actor);
```

---

## Vision

I have $$\gamma\colon [t_0, t_f] \to \mathbb{R}^{3}$$. I can render a space curve via

```cpp
auto actor1 = vtkm::tube_curve(gamma, t0, tf);
auto actor2 = vtkm::frenet_frame(gamma, t0, tf);
auto image = vtkm::render(actor1, actor2, ...);
```

![left](figures/FrenetFrame.png)

---

References:

- [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html)

- Moreland, Kenneth. "Diverging color maps for scientific visualization." International Symposium on Visual Computing. Springer, Berlin, Heidelberg, 2009.

- Moreland, Kenneth. "Why we use bad color maps and what you can do about it." Electronic Imaging 2016.16 (2016): 1-6.

---

- Pharr, Matt, Wenzel Jakob, and Greg Humphreys. "Physically based rendering: From theory to implementation." Morgan Kaufmann, 2016.

- Barr, Alan H. "Ray tracing deformed surfaces." ACM SIGGRAPH Computer Graphics 20.4 (1986): 287-296.

- Palais, Richard S. "The visualization of mathematics: towards a mathematical exploratorium." Notices of the AMS 46.6 (1999): 647-658.

---

- Pressley, Andrew N. "Elementary differential geometry." Springer Science & Business Media, 2010.

- Press, William H., et al. "Numerical recipes in C++." The art of scientific computing 2 (1992): 1002.

- Scott, Melvin, Beny Neta, and Changbum Chun. "Basin attractors for various methods." Applied Mathematics and Computation 218.6 (2011): 2584-2599.

---

- Yalcin, Hulya, Mustafa Unel, and William Wolovich. "Implicitization of parametric curves by matrix annihilation." International Journal of Computer Vision 54.1 (2003): 105-115.

- Cuyt, Annie AM, and Louis B. Rall. "Computational implementation of the multivariate Halley method for solving nonlinear systems of equations." ACM Transactions on Mathematical Software (TOMS) 11.1 (1985): 20-36.

- Barr, Alan H. "Ray tracing deformed surfaces." ACM SIGGRAPH Computer Graphics 20.4 (1986): 287-296.

---

- Park, Taezoon, Joonghyun Ji, and Kwang Hee Ko. "A second order geometric method for ray/parametric surface intersection." Computer Aided Geometric Design 30.8 (2013): 795-804.

- Axel Huebl, David Pugmire, Felix Schmitt, Richard Pausch, and Michael Bussmann. “Visualizing the Radiation of the Kelvin–Helmholtz Instability.” IEEE Transactions on Plasma Science 42 (10), October 2014.

- Toth, Daniel L. "On ray tracing parametric surfaces." ACM SIGGRAPH Computer Graphics 19.3 (1985): 171-179.