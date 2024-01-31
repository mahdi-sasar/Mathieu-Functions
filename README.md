# Update
The code presented here has now been merged into [BBN-Q's Mathieu function library](https://github.com/BBN-Q/MathieuFunctions.jl) and is accessible via Julia lang's general registry.
# Mathieu-Functions
This is a Julia language library for computing periodic Mathieu functions with arbitrary q
and n values. Mathieu functions are a class of very interesting and useful, but incredibly
finnicky equations. They come about as the solutions to the following differential
equation:

$$\frac{d^2}{dx^{2}} \psi(x) +(a-2q\cos(2x))\psi(x)=0$$

The specific problem that this package deals with is as follows: I need to plot and use
**periodic** Mathieu functions with unusually large $q$ values. Unusually large here means
$q$ values for which an advanced mathematics software like Mathematica, which is what
comes in handy in doing high level formal mathematics in DFT, quantum mechanics and
physics, **cannot** plot this function.

The word **periodic** is also key here, as Mathieu equation admits non-periodic solutions
that I am not interested in at this point for my study. I might add them here as well in
the future, but the periodic eigenfunctions are quite enough at this stage.

In most cases, we deal with small $q$ values and calculations with whatever computational
software you use can proceed quite smoothly. But for a particular research paper, I need
these functions for large negative and positive $q$ values, which comes about in the
analysis of a real life quantum system.

(If the related paper is published, I will update the read me files to include a reference
to it as well, to those who are interested.)

I have rigged a simple Julia package to calculate characteristic functions and periodic
Mathieu functions. It is blazingly fast and therefore, it allows expansion of other
periodic functions in terms of these by allowing fast numerical integrations of these
functions. Something that is very tricky again for large negative or positive $q$ values
in Mathematica. It can, of course, be used for all values of the $q$ parameter with the
added benefit of beating Mathematica in plotting them for $q$ values that stymies
Mathematica.

After writing these functions, I realized that there are some other Julia libraries
already available that calculate the characteristic values for difference $q$ values as
well. Please be sure to check them out as well: 
* [BBN-Q's Mathieu function library](https://github.com/BBN-Q/MathieuFunctions.jl)
* [jebej's Mathieu function library ](https://github.com/jebej/Mathieu.jl)

My calculation follows the brilliant textbook of Morse & Feshbach and the [*Digital
Library of Mathematical Functions*](https://dlmf.nist.gov/28.2). The great efficiency of
the Julia language gave rise to a pretty handy tool here that helped me greatly in my own
calculations. I hope it can aid you in your projects as well!


As time goes on, I will pour in more technical details about the project and how the whole
thing was developed in the end. It is a work in progress. So be sure to check in from time
to time and leave a star on the repo, if you have found this useful.


