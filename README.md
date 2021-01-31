slidenumbers: true

## Intro to Raytracing Differentiable Objects

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


![inline](figures/gaussian_curvature.png)

---

References:

- [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html)

- Moreland, Kenneth. "Diverging color maps for scientific visualization." International Symposium on Visual Computing. Springer, Berlin, Heidelberg, 2009.

- Moreland, Kenneth. "Why we use bad color maps and what you can do about it." Electronic Imaging 2016.16 (2016): 1-6.