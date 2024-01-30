# Mathieu-Functions
This is a Julia language library for computing periodic Mathieu functions with arbitrary q
and n values. Mathieu functions are a class of very interesting and useful, but incredibly
finnicky equations. They come about as the solutions to the following differential
equation:
$$
\begin{equation}
\frac{d^2}{dx^{2}} \psi(x) +(a-2q\cos(2x))\psi(x)=0
\end{equation}
$$
The problem that this package deals with is specifically as follows:
I need to plot and use periodic Mathieu functions with unusually large $q$
values. Unusually large $q$ values here means that advanced math softwares like
Mathematica, which is what comes in handy in doing high level formal mathematics in DFT,
quantum mechanics and physics, **cannot** plot this function.

In most cases, we deal with small $q$ values and the calculation can proceed quite
smoothly. But for a particular research paper, I need these functions for large $q$
values, which comes about in the analysis of a real life quantum system.

(If the paper is published, I will update the read me files to include a reference to it
as well, to those who are interested)

So what's wrong here? Why do softwares struggle with such plots?
The reason is as follows:
for periodic Mathieu functions of large $q$, the dependence of the solution to the exact
value of the characteristic value $a$ grows drastically. In fact, for $q=-60$, you need an
accurate $a$ with at least 10 decimal points. And this increases with increasing $q$.

I have rigged a very quick Julia package to calculate characteristic functions and
periodic Mathieu functions that does precisely that. It can, f course, be used for all
manner of $q$ parameters with the added benefit of beating Mathematica in plotting them
for $q$ values that stymies Mathematica.

This calculation follows the brilliant textbook of Morse & Feshbach (I mean, who
else?!). OMrse & Feshbach coupled with the great efficiency of the Julia language have
given rise to a pretty handy tool here that helped me greatly in my own calculations. I
hope it may aid you in yours!

**Problem** : Floats tend to under or over flow. Big arithmetic problem. The FiniteFloat
package seems to help out with this. It helps reduce the occurrences of NaNs and Infs,
which is exactly what we're looking for. This is a big problem for higher frequency
Mathieu functions.

**Main Problem**: Main issue is that for large $m$ and large $q$, the recursive formulae
are extremely unstable and tend to blow up. You'll need to *renormalize* the
characteristic function so that the recursive formulae **converge**.

# A Great Jump Ahead
The recursive formula for obtaining the coefficients of the Fourier expansion of the
Mathieu equations seems to blow up uncontrollably for high q values and even higher
indices.
I was able to by-pass this issue by using Legendre polynomials as a basis for expansion of
the Mathieu functions.
However, this has not been documented before and you are in a novel territory here.

