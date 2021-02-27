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

## History of raytracing

- Credited to Albrecht Dürer, but general purpose computing revitalized its study.

---

## Raytracing from a 1000 feet up

![inline](figures/RaysViewportSchema.png)

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

# Goal for Rendering Simulations

Render that data in the same way it's represented in the solver.

Barring that, render the data to the same order of accuracy that the solver produces.

---

## How to build a raytracer

- Read [Raytracing in a Weekend](https://raytracing.github.io/) by Peter Shirley.

- Compute intersection of ray $$\mathbf{P}(t) = \mathbf{O} + t\mathrm{D}$$ with parametric surface $$f(u,v)$$, or implicit surface $$f(x,y,z) = 0$$.

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

But let's examine a hack that makes it a bit better.

---

## "Patched" Newton's method: NR§9.7

Given $$\mathbf{f}\colon \Omega \subset \mathbb{R}^3 \to \mathbb{R}^3$$ defined as $$\mathbf{f}(t,u,v) = \mathbf{o} + t\mathbf{d} - \sigma(u,v)$$, we wish to find $$(t,u,v)$$ s.t. $$\mathbf{f}(t,u,v) = \mathbf{0}$$.

Newton update: $$(t_{k+1}, u_{k+1}, v_{k+1}) = (t_{k}, u_{k}, v_{k}) -\mathbf{J}_{\mathbf{f}}^{-1}\mathbf{f}(t,u,v)$$.

Patch: Define $$f(t,u,v) := \frac{1}{2}\left<\mathbf{f}(t,u,v), \mathbf{f}(t,u,v) \right>$$.

---

## Patched Newton's method

If $$f(t_{k+1}, u_{k+1}, v_{k+1}) > f(t_{k}, u_{k}, v_{k})$$ or $$(t_{k+1}, u_{k+1}, v_{k+1}) \not \in \Omega$$, reject the Newton step.

Define $$g(\lambda) := f(t_{k} + \lambda \delta t, u_{k} + \lambda \delta u, v_{k} + \lambda \delta v)$$, $$\lambda \in (0,1)$$.

Choose $$\lambda$$ that minimizes $$g$$.

---

![left](figures/helicoid_newton_wbacktracking.png)

Helicoid rendered with Newton's method patched with backtracking. Initial guess for $$t$$ provided by bounding cylinder intersection, no bounced rays. Note: Non-existence + Newton's method = Performance bug!

---

## Ray-tracers are 3D radiation transport codes!

What information do physically based renderers like [PBRT](http://www.pbr-book.org/3ed-2018/Shapes/Spheres.html) require for the intersection of a surface $$f(u,v)$$ with a ray?

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

## Implicit surfaces

Given a surface $$\mathcal{S} := \{ \mathbf{r} \colon f(\mathbf{r}) = 0\}$$, how do we generate ray intersections?

Naive idea: Use Newton's method:

$$t_{k+1} = t_{k} - \frac{\nabla f(r(t_k)) \cdot \mathbf{d}}{f(r(t_k))}$$

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
