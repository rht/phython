#0. List of papers?
 0.1 Studies of Nonlinear Problems
    1. what is the motivation of the paper? Why did Fermi bother to choose this
       system as the first ever numerical experiment?
    => http://arxiv.org/pdf/nlin/0411062v3.pdf "For Fermi this study was directly related to one his ﬁrst papers of 1923
       [2] in which he tried to rigorously prove the ergodicity hypothesis which
       lies at the core of traditional statistical mechanics"
    2. What is the progress of the paper?
    -> http://arxiv.org/pdf/nlin/0411062v3.pdf "About ten years later, two alternative explanations of the FPU paradox were suggested, giving rise to new phenomena, the integrability of nonlinear equations and dynamical (deterministic) chaos."
 0.2 The Fermi-Pasta-Ulam problem: Paradox turns discovery
 http://www.sciencedirect.com/science/article/pii/037015739290116H
 0.3 The story of Mary Tsingou
 http://arxiv.org/abs/0801.1590
 0.4 Nicest review: read this first before anything
 http://www.scholarpedia.org/article/Fermi-Pasta-Ulam_nonlinear_lattice_oscillations
 0.5 Nicest mathematical explanation
 http://www.scholarpedia.org/article/Fermi_Pasta_Ulam_systems_(FPU):_mathematical_aspects
 "In other words: the FPU chain does thermalize, but only very slowly."

(wordless:
Very few believed [localization] at the time, and even fewer saw its importance; among those who failed to fully understand it at first was certainly its author. It has yet to receive adequate mathematical treatment, and one has to resort to the indignity of numerical simulations to settle even the simplest questions about it.
--Philip W. Anderson, Nobel lecture, 8 December 1977
)

#Open questions
1. An important open question is what remains of the FPU phenomenology in two-dimensional and three-dimensional lattices.
   Worth to mention is the work of Benettin (Benettin 2004), who argues against the two-stage relaxation picture.
   http://chaos.aip.org/resource/1/chaoeh/v15/i1/p015108_s1?view=fulltext
   computing specific heat and thermal conductivity

#Questions
0. "The ergodic (german: ergodisch) behavior of such systems was studied with the primary aim of establishing, experimentally, the rate of approach to the equipartition of energy among the various degrees of freedom of the system"
-> What is ergodic behavior?
=> This is what drives unthermalised system with low entropy to a system with uniform temperature.

1. "this analysis amounts to a Lagrangian change of variables"
-> what does `Lagrangian` mean here?
=> It is used to exploit the fact that Lagrangian mechanism is coordinate-agnostic "all coordinates are created equal", while in Newton, if it is not a transformation from cartesian to another cartesian, then don't ever call it a "coordinate". so it's a Lagrangian-style change of coordinate.

2. "Instead of a gradual, continuous flow of energy from the first mode to the higher modes, all of the problems show an entirely different behavior."
-> I have seen this `physical prediction` before, but why do we expect this behavior in the first place?
=> let's see, why does the energy prefer to flow upward?
If you see the way they draw the sphere of density of states in 8.044: there are much more storage of states for higher eigenmodes.
What does this imply? Ultraviolet catastrophe, right?
=> (optional, recommended only if you don't mind wasting time...) ultraviolet catastrophe: you must have seen a lot about rayleigh-jeans formula being mentioned, but never derived. Hence the question "how do the people in the past derive the complicated formula of blackbody radiation FROM SCRATCH? without all the new ideas of QM?" it's in FloP ch 41 (my ground state of ideas is usually FloP).
=> for density of states in 1D in FloP, vol 3, page 4-10 and so on
"This result comes up again and again in many problems and should be memorized"
Yes, it is THAT IMPORTANT, equation 4.38
2.1
    -> you may wonder "many problems? like what?"
    => anything related to wave: in atomic physics, when dealing with the transition rate one particular molecular configuration into another, obviously in QED, in any random condensed matter system when calculating any density-of-state-related properties.

3. "The computation in the a_k variables  could have been more instructive for the purpose of observing directly the interaction between the a_k's."
-> Why do we bother to switch to a_k?
=> Because in a_k basis, the a_k's don't mix each other in the absence of nonlinearity


4. What is the cause behind the small wriggling oscillation?

5. Will you be able to observe thermalization?
scholarpedia: "It has been observed that relaxation to equipartition for long wavelength initial conditions proceeds in two stages. On a short time scale a packet of low frequency normal modes is formed, with the higher modes cut-off exponentially (Fucito et al. 1982): this phase is associated with the formation of the soliton train of Figure 3. On longer times, energy is slowly flowing to high frequency modes. This scenario has been recently advocated by Galgani and coworkers (Berchialla et al 2004) on the basis of numerical simulations and of a resonant normal form approach to the KdV equation (Bambusi and Ponno 2006)."
=> Yes! see the 300000 cycle run on sep 8. takes 54 min


6. Implement a travelling soliton!
=> you have to implement periodic boundary condition
   sep8, turns out we can't have a travelling soliton from fermi pasta ulam model yet. why?
   attempting to implement a vectorized version of the equation of motion


7. Can we calculate the Poincare cycle?

8. Quantify the recurrence period as a function of number of degrees of freedom.

9. What is the relation between kdV and FPU equation?
kdV is the continuum limit of FPU http://en.wikipedia.org/wiki/Korteweg%E2%80%93de_Vries_equation

#progress
#sep8
* found the original paper on soliton! from Zabusky's story
    http://chaos.aip.org/resource/1/chaoeh/v15/i1/p015102_s1?view=fulltext
    Interaction of "Solitons" in a Collisionless Plasma and the Recurrence of Initial States
    http://prl.aps.org/abstract/PRL/v15/i6/p240_1

#scicomp features
parallel computing library comparison:
1. parallel python: lots of dependencies
2. pprocess
http://www.astrobetter.com/parallel-processing-in-python/
pprocess vs parallel python vs multiprocessing vs mincemeat
