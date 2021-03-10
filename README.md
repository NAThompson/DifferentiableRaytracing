slidenumbers: true

## Raytracing Differentiable Surfaces

![inline](figures/ornl_logo.svg)

Nick Thompson

---

## What is raytracing

- Cast lights from eyes to objects, and see if they reflect into a light source
- Exactly the opposite of real world, but no problem due to electromagnetic reciprocity

---

![](figures/albrecht_durer.jpg)

---

## How to render a sphere

Compute the intersection of a line and a sphere. Easy.

Line: $$\mathbf{P}(t) = \mathbf{O} + t\mathbf{D}$$
Sphere: $$(\mathbf{P} - \mathbf{C})^2 = R^2$$.

Plug in to get quadratic equation for intersection points.

But isn't a sphere differentiable?

---

## Raytracing *differentiable* objects

Naive raytracing is *agnostic* to the differentiability of the objects in your scene.

All that matters to raytracing is that you can compute the intersection of your objects with a line.

---

## Ok, why do we need to focus on raytracing *differentiable* structures?

Raytracing triangle meshes has become dominant in game programming. The assumption is that *creation* of the mesh is not that expensive-a good assumption for game asset generation; a bad one for visualizing ODE/PDE solutions.

CAD programs *used* to render diagrams using smoother surfaces; e.g., Bezier and Catmull-Clark.

---

## A billion triangles a frame!

![inline](https://youtu.be/qC5KtatMcUw?t=106)

---

## Interpolation in the scalar case

Tri-meshes are a form of linear interpolation on $$\mathbb{R}^3$$.

If we restrict to 1D, what benefit is there to smoother representation?

---

The error of linear interpolation is $$\mathcal{O}(h^2)$$.

The error of a cubic B-spline is $$\mathcal{O}(h^4)$$.

So for the same level of error, a cubic B-spline requires a square root of the amount of memory that linear interpolation requires!

This requires differentiability of the data; if it's not differentiable, just use linear interpolation.

---


Again, rendering massive data is *not* a problem. It's *creating* that data in a simulation.

Try solving your ODEs with Euler's method, or PDEs with first-order accurate finite differences.

---

## Goal (dream?) for rendering simulations

Render the data in the same way it's represented in the solver.

If you do your simulation using Chebyshev spectral elements, VTK-m should be rendering your data using your Clenshaw recurrence.


---

## How to build a raytracer

- Read [Raytracing in a Weekend](https://raytracing.github.io/) by Peter Shirley.

- Compute intersection of ray $$\mathbf{P}(t) = \mathbf{O} + t\mathrm{D}$$ with parametric surface $$\sigma(u,v)$$, or implicit surface $$f(x,y,z) = 0$$.

- Determine axis aligned bounding box of objects to place it in the bounding volume hierarchy.

---

## Sphere -> Ellipsoid

Ellipsoids are just affine transformations of spheres. For an ellipsoid
$$\sum_{i = 0}^{2} \left( \frac{x_i - c_i}{s_i} \right)^2 = 1$$
and a ray $$P(t) = \mathbf{O} + t\mathrm{D}$$, the ray intersects at
$$
\sum_{i=0}^{2} \left( \frac{o_{i} + td_{i} - c_i}{s_i}  \right)^2 = 1
$$

Again, just a quadratic.

---

![inline](figures/ellipsoid.png)


---

Ok, but: In SciViz, we need to display scalar fields on 3D objects to represent properties or phenomenom we have no sensory apparatus for. This is called pseudocoloring.

e.g., flow fields, radiation, curvature.

So we need a color map: $$f\colon \mathbb{R} \to [0,255]\times [0,255] \times [0,255]$$.

---

In movies and game engines, colors are chosen to maximize aesthetics.

There's nothing wrong with aesthetics in sciviz, but our first allegiance is to scientific insight.

---

## Criteria for SciViz colormaps

> The efficacy of a pseudocolor visualization is contingent on the ability of a human observer to translate the colors back into the numeric values they represent.
-- Ken Moreland

---

## Sensible color maps for raytraced surfaces.

See [kennethmoreland.com/color-advice](https://www.kennethmoreland.com/color-advice/)

![inline](figures/plasma_map.png)

---

## Ellipsoid colored by Gaussian Curvature

![inline](figures/gaussian_curvature.png)

---

## Parametric Equation of a Torus

$$(x^2 + y^2 + z^2 + R^2 - r^2)^2 - 4R^2(x^2+y^2) = 0$$

Ray $$\mathbf{O} + t\mathbf{D}$$ intersects torus when $$x = o_x + td_x, y = o_y + td_y, z = o_z + td_z$$, which gives a quartic equation for $$t$$.

- No real roots implies no intersection of ray with torus
- Multiple real roots implies multiple intersections; take the smallest real root (point of surface closest to camera)

---

## Gaussian Curvature of Torus

![inline](figures/torus_curvature.png)


---

## Gaussian Curvature of Torus in Paraview

![inline](figures/pv_torus_curvature.png)

---

## Parametric equation of helicoid

$$\sigma(u, v) = (rv\cos(2\pi u), rv\sin(2\pi u), \lambda (u-1/2))$$.

Intersects ray when

$$
o_x + td_x = rv\cos(2\pi u),
o_y + td_y = rv\sin(2\pi u),
o_z + td_z = \lambda(u-1/2)
$$

---

## Helicoid ray intersections

Newton's method will work well enough here, but getting sensible starting values is key.

Bound the helicoid with a cylinder. A cylinder is a quadric surface so we can obtain $$[t_{\min}, t_{\max}]$$ by a quadratic equation.

Then $$r^2v^2 = (o_x + td_x)^2 + (o_y + td_y)^2$$ gives us bounds on $$[v_{\min}, v_{\max}]$$, and $$u = (o_z+td_z)/\lambda + 1/2$$ gives us bounds on $$[u_{\min}, u_{\max}]$$.

---

## Inequalities are just as important as numerics!

Working to *constrain* a ray intersection with an object is just as important as the convergence of the method to the correct solution.

The ability to say "this ray *does not* intersect our object" in as few flops as possible is just as important for your performance as the numerical method to compute the intersection when it does.

---

The Newton fractal teaches us that it's difficult to predict *which* root Newton's method with a given starting point will converge to.

In ray tracing, we need the *minimal* root to determine the closest intersection point to the ray.

![left](figures/newton_fractal.jpg)

^ A Newton fractal colored by root reached.

---



![left](figures/wrong_solutions.png)

A helicoid pseudocolored by Gaussian curvature. The pseudocoloring changes discontinuously where the non-minial roots are identified.

---

## Multivariate root-finding is the worst!

![inline](figures/the_worst.gif)


---

## Multivariate root-finding is the worst!

> We make an extreme, but wholly defensible, statement: There are *no* good, general methods for solving systems of more than one nonlinear equation.
-- Numerical Recipes

---

## Multivariate root-finding is the worst!

Bounced rays make the problem even worse than the generic multivariate rootfinding problem.

A ray that bounces off an object satisfies the equation $$f(\mathbf{o} + t\mathbf{d}) = 0$$ at $$t = t_{\min}$$, and it *might* hit the object again.

So you need to find the $$t$$ *closest* to $$t_{\min}$$ that hits the object, but not the $$t = t_{\min}$$ solution.

---

And numerical imprecision might mean your $$t = t_{\min} + \epsilon$$ solution is essentially the *same* as your $$t = t_{\min}$$ solution.

This generates a visual artefact known as [shadow acne](https://digitalrune.github.io/DigitalRune-Documentation/html/3f4d959e-9c98-4a97-8d85-7a73c26145d7.htm).

![left](figures/shadow_acne.jpg)

---

![inline](figures/bad_mkay.jpg)

---

## Newton's method is bad

- Need to compute derivatives/Jacobians
- Sensitive to rounding errors
- Difficulty with multiple roots
- There can exist full-measure regions of initial conditions where it diverges.
- The Newton iterate can converge, but to a value which is not a root.
- Converges to the "wrong root"

But let's examine a hack that makes it a bit better.

---

## "Patched" Newton's method: NR§9.7

Given $$\mathbf{f}\colon \Omega \subset \mathbb{R}^3 \to \mathbb{R}^3$$ defined as $$\mathbf{f}(u,v,t) = \sigma(u,v) - (\mathbf{o} + t\mathbf{d})$$, we wish to find $$(u,v, t)$$ s.t. $$\mathbf{f}(u,v,t) = \mathbf{0}$$.

Newton update: $$(u_{k+1}, v_{k+1}, t_{k+1}) = (u_{k}, v_{k}, t_{k}) -\mathbf{J}_{\mathbf{f}}^{-1}\mathbf{f}(u,v,t)$$.

Patch: Define $$f(u,v,t) := \frac{1}{2}\left<\mathbf{f}(u,v,t), \mathbf{f}(u,v,t) \right>$$.

---

## Patched Newton's method

If $$f(u_{k+1}, v_{k+1}, t_{k+1}) > f(u_{k}, v_{k}, t_k)$$ or $$(u_{k+1}, v_{k+1}, t_{k+1}) \not \in \Omega$$, reject the Newton step.

Define $$g(\lambda) := f(u_{k} + \lambda \delta u, v_{k} + \lambda \delta v, t_{k} + \lambda \delta t)$$, $$\lambda \in (0,1)$$.

Choose $$\lambda$$ that minimizes $$g$$.

---

![left](figures/helicoid_newton_wbacktracking.png)

Helicoid rendered with Newton's method patched with backtracking. Initial guess for $$t$$ provided by bounding cylinder intersection, no bounced rays.

Note: Non-existence + Newton's method = Performance bug!

---

## Performance

Since we are restricted to $$n \le 4$$ in computer graphics, we can soberly consider rootfinding methods that are asymptotically ill-advised.

e.g., multivariate Halley's method, which requires formation of Hessian tensors.

---

## Halley's method

$$t_{k+1} = t_{k} - \frac{2f(t_k)f'(t_k)}{2f'(t_k)^2 - f(t_k)f''(t_k)}$$

Convergence is cubic, much more importantly, *which* root Halley's method converges to is much more predictable; see [here](http://faculty.nps.edu/bneta/papers/BasinsThirdFourth.pdf) and [here](http://dx.doi.org/10.1016/j.amc.2011.07.076).

---

![left](figures/z5m1_newton.jpg)
![right](figures/z5m1_halley.jpg)

^ [Newton vs Halley fractals](http://www.hpdz.net/StillImages/Newton-Halley.htm)

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

The $$\otimes$$ and $$\oslash$$ are bizarre componentwise multiplications and divisions; see [Cuyt](https://dl.acm.org/doi/pdf/10.1145/3147.3162) for details.


---

## Implicitization

Ray intersection with a parametric surface is a multivariate rootfinding problem.

Ray intersection with an implicit surface is a scalar rootfinding problem-much easier! Non-existence characterization simpler, fallback to bisection possible.

Consider [converting](https://www.cs.cmu.edu/~hulya/Publications/IJCV03Paper.pdf) your parametric surfaces into implicit surfaces!

---

## So many more techniques

- Use hit points from adjacent pixels as starting values for Newton/Halley iterates.
- Using fundamental forms, change into a coordinate system where rays are bent and surfaces are flat; see [Barr](https://dl.acm.org/doi/pdf/10.1145/15886.15918).
- Do a coarse triangulation of the surface, then use a Halley or Newton iterate to refine the hit points off the triangulation.
- Bounding volume hierarchies.

---

## So many more techniques

- Is $$\mathbf{f}(u,v,t) = \sigma(u,v) - \mathbf{o} - t \mathbf{d}$$ contractive? Banach fixed point theorem.
- Load balancing and parallelization: Discussed in [Physically Based Rendering](http://www.pbr-book.org/3ed-2018/Utilities/Parallelism.html).


---

## Rendering space curves

The probability that a space curve $$\gamma \colon [t_0, t_f] \to \mathbb{R}^{3}$$ intersects a ray $$\mathbf{o} + t\mathbf{d}$$ is zero!

Space curves must be extended, either via glyphs, or as "tubes". Constructing a tube from a space curve $$\gamma$$:

$$
\sigma(\theta, t) = \gamma(t) + r(\mathbf{n}(t)\cos(\theta) + \mathbf{b}(t) \sin(\theta))
$$

where $$\mathbf{n}$$ is the principle normal, and $$\mathbf{b}$$ is the binormal.

---

## Ray-tracers are 3D radiation transport codes!

What information do physically based renderers like [PBRT](http://www.pbr-book.org/3ed-2018/Shapes/Spheres.html) require for the intersection of a surface $$\sigma(u,v)$$ with a ray?

- The intersection point $$\mathbf{p}$$, as well as $$t$$ and $$(u,v)$$.
- The normal to the surface $$\mathbf{n}$$.
- Partial derivatives $$\partial_{u}\mathbf{p}$$, $$\partial_{v}\mathbf{p}$$, $$\partial_{u}^2\mathbf{p}$$, $$\partial_{v}^2\mathbf{p}$$, $$\partial_{uv}^2\mathbf{p}$$
- The first and second [fundamental forms](https://en.wikipedia.org/wiki/First_fundamental_form).

---

## Ray-tracers are 3D radiation transport codes!

- The information required for radiation transport is a pretty good base for doing generic sciviz!

- Differential geometry is a solid foundation to build your SciViz on.

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

I have $$\gamma\colon [t_0, t_f] \subset \mathbb{R} \to \mathbb{R}^{3}$$. I can render a space curve via

```cpp
auto actor1 = vtkm::tube_curve(gamma, t0, tf);
auto actor2 = vtkm::frenet_frame(gamma, t0, tf);
auto image = vtkm::render(actor1, actor2, ...);
```

![left](figures/FrenetFrame.png)

---

## Is there a path from meshes to smooth surfaces?

Yes, [subdivision surfaces](http://www.pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces.html).

These form smooth surfaces out of mesh points.

---

References:

- [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html)

- Moreland, Kenneth. "Diverging color maps for scientific visualization." International Symposium on Visual Computing. Springer, Berlin, Heidelberg, 2009.

- Moreland, Kenneth. "Why we use bad color maps and what you can do about it." Electronic Imaging 2016.16 (2016): 1-6.

- Pharr, Matt, Wenzel Jakob, and Greg Humphreys. "Physically based rendering: From theory to implementation." Morgan Kaufmann, 2016.

---

- Barr, Alan H. "Ray tracing deformed surfaces." ACM SIGGRAPH Computer Graphics 20.4 (1986): 287-296.

- Palais, Richard S. "The visualization of mathematics: towards a mathematical exploratorium." Notices of the AMS 46.6 (1999): 647-658.

- Pressley, Andrew N. "Elementary differential geometry." Springer Science & Business Media, 2010.

- Press, William H., et al. "Numerical recipes in C++." The art of scientific computing 2 (1992): 1002.

---

- Scott, Melvin, Beny Neta, and Changbum Chun. "Basin attractors for various methods." Applied Mathematics and Computation 218.6 (2011): 2584-2599.

- Yalcin, Hulya, Mustafa Unel, and William Wolovich. "Implicitization of parametric curves by matrix annihilation." International Journal of Computer Vision 54.1 (2003): 105-115.

- Cuyt, Annie AM, and Louis B. Rall. "Computational implementation of the multivariate Halley method for solving nonlinear systems of equations." ACM Transactions on Mathematical Software (TOMS) 11.1 (1985): 20-36.

- Barr, Alan H. "Ray tracing deformed surfaces." ACM SIGGRAPH Computer Graphics 20.4 (1986): 287-296.

---

- Park, Taezoon, Joonghyun Ji, and Kwang Hee Ko. "A second order geometric method for ray/parametric surface intersection." Computer Aided Geometric Design 30.8 (2013): 795-804.

- Axel Huebl, David Pugmire, Felix Schmitt, Richard Pausch, and Michael Bussmann. “Visualizing the Radiation of the Kelvin–Helmholtz Instability.” IEEE Transactions on Plasma Science 42 (10), October 2014.