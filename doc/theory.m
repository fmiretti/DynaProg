%[text] # Theory
%[text] Explore the theory behind the dynamic programming algorithm implemented in DynaProg
%%
%[text] ## Introduction
%[text] DynaProg solves finite horizon, deterministic multi-stage decision problems.
%[text] A multi-stage decision problem is a control problem where decisions must be made in stages in order to minimize a certain cost. A system, which is characterized by its state, evolves through the stages and this evolution is influenced by the decisions themselves. 
%[text] Deterministic problems in particular refer to the fact that the system's evolution for each stage is fully predictable: given the current state and the decision variables, the state of the system at the next stage can be evaluated with no uncertainty.
%[text] Finally, a finite-horizon problem is a control problem where decisions must be made for a finite number of stages.
%[text] DynaProg solves these problems using a dynamic programming algorithm. Dynamic programming was first introduced by Richard Bellmann \[1\], and has since been applied in a wide variety of contexts and its theoretical basis extended by many other great contributors.
%[text] ## **Basics**
%[text] Given the initial state of the system $x\_0${"altText":"x\_0"}, our control optimization problem consists in selecting the sequence of control variables over$N${"altText":"N"}stages $(u\_o, u\_1, ..., u\_{N-1})${"altText":"(u\_o, u\_1, ..., u\_(N-1))"} that minimizes some total cost 
%[text]  $J(x\_0, u\_0, u\_1, ..., u\_{N-1}) = F(x\_N) + \\sum\_{k=0}^{N-1} L\_k(x\_k,u\_k),$
%[text] where $ L\_k(x\_k,u\_k)${"altText":" L\_k(x\_k,u\_k)"} is the stage cost incurred in advancing by one stage and $F(x\_N)${"altText":"F(x\_N)"} is a terminal cost associated to the terminal system state.
%[text] In other words, the goal of the optimization problem is therefore to find the optimal cost $V\_0(x\_0)$, that is the minimum total cost that can be incurred:
%[text] $V\_0(x\_0) = \\min\_{\\matrix{u\_k \\in U\_k(x\_k),\\cr k=1,...,N-1}} \\left( F(x\_N) + \\sum\_{k=0}^{N-1} L\_k(x\_k,u\_k) \\right),$
%[text] and the control sequence that minimizes it:
%[text] ${u^\*\_0,...,u^\*\_{N-1}} = \\underset{\\matrix{u\_k \\in U\_k(x\_k),\\cr k=1,...,N-1}}{\\mathrm{argmin}} \\left( F(x\_N) + \\sum\_{k=0}^{N-1} L\_k(x\_k,u\_k) \\right).$
%[text] Here, $U\_k(x\_k)${"altText":"U\_k(x\_k)"} is the set of feasible control variables that can be selected at stage $k${"altText":"k"}.
%[text] The foundation of dynamic programming-based optimization algorithms is Bellman's principle of optimality:
%[text] *"An optimal policy has the property that whatever the initial state and  initial decision are, the remaining decisions must constitute an optimal policy with regard to the state resulting from the first decision"* \[1\].
%[text] In practice, the idea is to divide the original control problem into several tail subproblems, each of which is defined as the optimal control problem starting from a certain stage $k${"altText":"k"} to the final stage $N${"altText":"N"}, and to iteratively solve increasingly long tail subproblems by exploiting the solution of the previously obtained tail subproblems.  
%[text] The optimization is divided in two phases: a so-called "backward" phase and a "forward" phase. 
%[text] ### DP algorithm (value-based)
%[text] #### Backward Phase
%[text] Start by setting:
%[text] $V\_N(x\_N) = F(x\_N).$
%[text] Set $k = N-1$and solve:
%[text] $V\_k(x\_k) = \\min\_{u\_k \\in U\_k(x\_k)} \\left( L\_k(x\_k,u\_k) + V\_{k+1}\\left( f\_k(x\_k,u\_k)\\right) \\right) \\quad (1)$
%[text] for all $x\_k${"altText":"x\_k"}.
%[text] Set $k$ back by one stage and repeat the last step until $V\_0(x\_0)$ is obtained. 
%[text] Because $V\_k(x\_k)$ represents the minimum cost that can be incurred if the system must evolve from stage $k${"altText":"k"} to stage $N${"altText":"N"} given the initial state $x\_k${"altText":"x\_k"}, it is also called the cost-to-go. Some sources also refer to this quantity as the value function.
%[text] #### Forward Phase
%[text] Start from $k=0${"altText":"k=0"} with an initial state $x\_0$, and evaluate
%[text] $u^\*\_k(x\_k) =  \\underset{u\_k \\in U\_k(x\_k)}{\\mathrm{argmin}} \\left( L\_k(x\_k,u\_k) + V\_{k+1}\\left( f\_k(x\_k,u\_k)\\right) \\right)$.
%[text] Advance the simulation by updating the state variables
%[text] $x\_{k+1} = f(x\_k, u^\*\_k(x\_k))${"altText":"x\_(k+1) = f(x\_k, u^\*\_k(x\_k))"}.
%[text] Advance $k${"altText":"k"} by one stage and repeat these two steps until the last stage.
%[text] ### Policy-based algorithm
%[text] By default, DynaProg uses a value-based algorithm, as described previously. Alternatively, a policy-based algorihm can be used. 
%[text] In this context, a policy is defined as a sequence of functions $\\mu\_k(x\_k)$ that directly gives the optimal control for stage $k$ if the state is equal to $x\_k${"altText":"x\_k"}.
%[text] In the backward phase, rather than storing the sequence of value functions$V\_k(x\_k)$, we directly store the controls that minimize $(1)$. With these values, we construct the function $\\mu\_k(x\_k)$ as a linear interpolant. Then, in the forward phase, we simply compute the optimal controls using the policy that was generated by the backward phase.
%[text] #### Backward Phase
%[text] Start by setting:
%[text] $V\_N(x\_N) = F(x\_N).$
%[text] Set $k = N-1$and solve:
%[text] $u^\*\_k(x\_k) =  \\underset{u\_k \\in U\_k(x\_k)}{\\mathrm{argmin}}  \\left( L\_k(x\_k,u\_k) + V\_{k+1}\\left( f\_k(x\_k,u\_k)\\right) \\right)$
%[text] for all $x\_k${"altText":"x\_k"}. Use these values to create a function $\\mu\_k(x\_k)$. 
%[text] Set $k$ back by one stage and repeat the last step until $\\mu\_0(x\_k)$ is obtained. 
%[text] #### Forward Phase
%[text] Start from $k=0${"altText":"k=0"} with an initial state $x\_0$, and evaluate
%[text] $u^\*\_k(x\_k) =  \\mu\_k(x\_k)$.
%[text] Advance the simulation by updating the state variables
%[text] $x\_{k+1} = f(x\_k, u^\*\_k(x\_k))${"altText":"x\_(k+1) = f(x\_k, u^\*\_k(x\_k))"}.
%[text] Advance $k${"altText":"k"} by one stage and repeat these two steps until the last stage.
%[text] ## **Discretization**
%[text] DynaProg applies the outlined algorithm to derive a numerical solution to the control optimization problem.  For any multi-stage process with discrete state and control variables, the solution is straightforward and it does not present any particular numerical hazard. If process instead is a continuous time process and/or the state and control variables are continuous, a numerical solution can still be achieved, but some care must be taken.
%[text] #### Time discretization
%[text] If the process considers a system which evolves in continuous time, then time must be discretized with a certain time step into a sequence of stages. An optimal solution can be then achieved for the discrete-time equivalent process.
%[text] #### Control variables discretization
%[text] If one or more of control variables are continuous, they must be discretized. This simplifies the numerical solution as the $\\mathrm{min}${"altText":"min"} and $\\mathrm{argmin}${"altText":"argmin"} operations required by the algorithm reduce to finding minimum values over finite sets. 
%[text] If the value-based algorithm is used, the optimal control sequence found by the optimization algorithm will be restricted to the discretized control variables. 
%[text] If the policy-based algorithm is used instead, the optimal control sequence evaluated in the forward phase will be evaluated by linear interpolation on the policy functions $\\mu\_k(x\_k)$. Nonetheless, the discrete control grids are used in the backward phase to generate these functions.
%[text] #### State variables discretization
%[text] If the system has one or more state variables, it does not need to be discretized in the simulation. 
%[text] However, the backward phase of the optimization algorithm relies on storing, at each iteration, $V\_k(x\_k)${"altText":"J^\*\_k(x\_k)"} (the value of the cost-to-go as a function of $x$). In order to do this, DynaProg evaluates $V\_k(x\_k)${"altText":"J^\*\_k(x\_k)"} for a finite set of sampled values of $x\_k${"altText":"x\_k"}: for this reason, a discretized computational grid for the state variables is needed. Then, when values of $V\_{k+1}\\left( f\_k(x\_k,u\_k)\\right) $ must be evaluated for values of $x\_{k+1} = f\_k(x\_k,u\_k)$ that do not belong to the computational grid, it is obtained by linear interpolation.
%[text] The discretization of this computational grid for the state variables affects the computation of the costs-to-go and, as such, it is a source of sub-optimality. Selecting the proper discretization level for the state variables grid usually requires some understanding of the physics behind the system under analysis and possibly some trial-and-error.
%[text] To wrap up, in the presence of contiunous state variables, a discrete computational grid must be created for the optimization algorithm. This has an influence on the optimality of the solution, even though the state is not really discretized in simulating the system's evolution. The trajectory of the state variables given the optimal control sequence is continuous and it is not tied to the discrete computational grid.
%[text] %[text:anchor:H_6C2F7696] ## **The terminal cost** 
%[text] DynaProg allows you to specify the terminal cost $F(x\_N)$ as a function handle by setting the `TerminalCost` property. Read the [Syntax](file:./syntax.mlx) guide for more information.
%[text] In addition to that, DynaProg adds a penalty term $\\Psi(x\_N)$ to the terminal cost in order to enforce terminal state constraints (if present). In other words, the cost-to-go at stage $N\n$ is initialized to:
%[text] $V\_N(x\_N) = F(x\_N) + \\Psi(x\_N)$.
%[text] Currently, two alternatives can be selected to define this penalty term: *rift* penalization and *linear* penalization.
%[text] %[text:anchor:M_BBA8F63A] *Rift* penalization means that the penalty term has the form:
%[text] $\\begin{cases}\n			\\Psi(x\_N) = V\_\\infty \\quad \\text{if}\\ x\_N \< x\_\\mathrm{lb} \\lor x\_N \> x\_\\mathrm{ub},\\\\\n			\\Psi(x\_N) = 0 \\quad \\text{otherwise.}\\\\\n	\\end{cases}$
%[text] where $ x\_\\mathrm{lb}$ and $x\_\\mathrm{ub}$ are the upper and lower bounds of the terminal state constraints. The term $V\_\\infty$ can be modified by the user using the `myInf` property.
%[text] *Linear* penalization means that the terminal cost has the form
%[text] $\\begin{cases}\n		\\Psi(x\_N) = V\_\\infty \\quad \\text{if }\\ x\_N \< x\_\\mathrm{lb} \\lor x\_N \> x\_\\mathrm{ub},\\\\\n		\\Psi(x\_N) = p\_\\Psi^T \\, \\min\\left( (x\_N - x\_\\mathrm{lb})^{| \\cdot |}, (x\_N - x\_\\mathrm{ub})^{| \\cdot |} \\right) \\quad \\text{otherwise,}\\\\\n\\end{cases}\n$
%[text] i.e. the terminal cost is proportional to the distance between the terminal state and the closest bound of the terminal state constraints, and $b$ is the proportionality factor. Here, $v^{| \\cdot |}$ denotes element-wise absolute value of all elements of the vector $v$.
%[text] The *linear* penalization can fail to enforce terminal state constraints if the $b$ factors are not set properly. *Rift* penalization on the other hand can cause the optimization to fail if the terminal state bounds are too tight or if the state grids are too coarse. This is caused by numerical issues related to the interpolations performed by the algorithm to obtain values for the cost-to-go function in the proximity of the unfeasible region.
%[text] You can specify what method to use by using the `VFPenalty` property and the proportionality factors $b$ by using the `VFFactors` property. Upper and lower bounds for the terminal state constraints are specified with the `StateFinal` property. Read the [Syntax](file:./syntax.mlx) guide for more information.
%[text] %[text:anchor:M_LevelSet] ## **The Level-Set method**
%[text] As an alternative to penalizing the terminal cost, DynaProg can enforce terminal state constraints with the *Level-Set* method \[2\]. You can enable it by setting the `UseLevelSet` property to `true`.
%[text] The idea is to decouple the enforcement of the terminal state constraints from the minimization of the cost. To this end, an auxiliary *level-set function* $I\_k(x\_k)$ is introduced alongside the cost-to-go $V\_k(x\_k)$ to keep track of which states can still reach the *target set*, i.e. the set of terminal states that satisfy the constraints, $\\mathcal{T} = \\lbrace x\_N : x\_\\mathrm{lb} \\le x\_N \\le x\_\\mathrm{ub} \\rbrace$.
%[text] At the final stage, the level-set function is initialized so that it is non-positive inside the target set and positive outside of it:
%[text] $I\_N(x\_N) = \\max\_i \\, \\max\\left( x\_{\\mathrm{lb},i} - x\_{N,i}, \\; x\_{N,i} - x\_{\\mathrm{ub},i} \\right),$
%[text] where the index $i$ spans the state variables. By definition, $I\_N(x\_N) \\le 0$ if and only if $x\_N \\in \\mathcal{T}$(meets the terminal constraints).
%[text] In the backward phase, the level-set function is propagated together with the cost-to-go:
%[text] $I\_k(x\_k) = \\min\_{u\_k \\in U\_k(x\_k)} I\_{k+1}\\left( f\_k(x\_k,u\_k) \\right).$
%[text] A state $x\_k$ at stage $k$ can reach the target set if and only if $I\_k(x\_k) \\le 0$. If $I\_k(x\_k) \>0$, that means that starting from state $x\_k$ there is no feasible sequence of controls $(u\_k, u\_{k+1}, ..., u\_{N-1})$ that can drive the state within the target set at stage $N$.
%[text] At each stage, the controls that lead to a reachable successor state define the *reachable control set* $U\_k^R(x\_k) = \\lbrace u\_k \\in U\_k(x\_k) : I\_{k+1}\\left( f\_k(x\_k,u\_k) \\right) \\le 0 \\rbrace \\subseteq U\_k(x\_k)$.
%[text] The cost-to-go update is restricted to the reachable control set, so that the minimization only considers controls that keep the terminal state constraints satisfiable:
%[text] $V\_k(x\_k) = \\min\_{u\_k \\in U\_k^R(x\_k)} \\left( L\_k(x\_k,u\_k) + V\_{k+1}\\left( f\_k(x\_k,u\_k) \\right) \\right).$
%[text] If the reachable control set $U\_k^R(x\_k)$ is empty, the target set cannot be reached from $x\_k$; in this case the cost-to-go is evaluated at the control that minimizes the level-set function, i.e. the one that comes closest to reaching the target set \[2\].
%[text] Unlike the value function penalization approches, the Level-Set method enforces the terminal state constraints without introducing a large penalty term near the unfeasible region. This makes it more robust against the numerical issues caused by interpolating the cost-to-go close to the constraint bounds.
%[text] ## Notes on computational time
%[text] Dynamic programming is an extremely powerful optimization tool in that it allow to solve problems even if the state dynamics and cost function are highly nonlinear in their structure. However, this comes at the cost of relatively high computational time. By far, the most time-consuming task is to building the cost-to-go function in the backward phase. Consider the cost-to-go update that is repeated for each stage:
%[text] $V\_k(x\_k) = \\min\_{u\_k \\in U\_k(x\_k)} \\left( L\_k(x\_k,u\_k) + V\_{k+1}\\left(  f\_k(x\_k,u\_k)\\right) \\right).$
%[text] Note that this requires evaluating $f\_k(x\_k,u\_k)${"altText":"f\_k(x\_k,u\_k)"}for each combination of state and control variables defined by their discrete computational grids. This means that, for each stage,  the state dynamics must be evaluated $N\_{SV\_1} \\times N\_{SV\_2} \\times ... \\times N\_{SV\_\_{NS}} \\times N\_{CV\_1} \\times N\_{CV\_2} \\times ... \\times N\_{CV\_{NC}}${"altText":"N\_{SV\_1} x N\_{SV\_2} x ... x N\_{SV\_\_{NS}} x N\_{CV\_1} x N\_{CV\_2} x ... x N\_{CV\_{NC}}"} times, where $N\_{SV\_i}${"altText":"N\_SV\_i"} and $N\_{CV\_j}${"altText":"N\_CV\_j"} are the number of grid points that discretize the state variable $SV\_i$ and the control variable $CV\_j${"altText":"CV\_j"} and$NS${"altText":"NS"} and $NC${"altText":"NC"}are the number of state and control variables.
%[text] Not only is the number of function evaluations that must be performed very high, but it also increases dramatically by adding more state and/or control variables. For this reason, for complex system dynamics and cost functions, it is important that no unneccesary computations are defined in the system and cost function by the user. DynaProg provides several ways to help reducing them to a minimum, with simple features such as [exogenous inputs](file:./handle_exogenous_inputs.mlx) and [additional inputs](file:./use_additional_inputs.mlx) and the more advanced [system and cost function split](file:./splitting_the_model.mlx).
%[text] ## References
%[text] \[1\] Richard Ernest Bellman, *Dynamic programming*. Princeton, Nj: Princeton University Press, 2010. ISBN: 9780691146683.
%[text] \[2\] P. Elbert, S. Ebbesen, and L. Guzzella, "Implementation of Dynamic Programming for n-Dimensional Optimal Control Problems With Final State Constraints," *IEEE Transactions on Control Systems Technology*, vol. 21, no. 3, pp. 924–931, 2013. doi: 10.1109/TCST.2012.2190935.
%%
%[text]{"align":"right"} *Contact: federico.miretti@polito.it*

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":34.1}
%---
