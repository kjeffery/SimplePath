# SimplePath

A simple path tracer for educational purposes. This path tracer has no ambition to become anything more than an
educational project. You are welcome to make changes and create pull requests.

[Physically Based Rendering](https://www.pbr-book.org/) is an amazing resource, but it may be too complicated to get
your feet wet. Many like [Ray Tracing in One Weekend](https://raytracing.github.io/), but I find some of the shortcuts
taken in the course of expediency to be bad habits and poor teaching methods. I think
[Ray Tracing from the Ground Up](https://www.oreilly.com/library/view/ray-tracing-from/9781498774703/) is a wonderful
teacher of concepts, but the code architecture leaves something to be desired.

Here I have tried to make concepts easy to follow, often at the cost of ray tracing performance or convergence rates. I
must admit that I have failed in some cases where I get too caught up in thinking about how things would be done in a
production path tracer. For instance, the vector classes use SIMD implementations, which is not necessary and makes the
code less portable. At the time, my thinking was that users could focus on the interface rather than the details.
