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

- Compute intersection of ray $$\mathbf{P}(t) = \mathbf{O} + t\mathrm{D}$$ with object: $$\mathbf{P}(t) = f(x,y,z)$$.

- Determine axis aligned bounding box of objects to place it in the bounding volume hierarchy.

---

References:

- [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html)