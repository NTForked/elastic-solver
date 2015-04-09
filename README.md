# FEM Solver #

This is an experimental implementation of the FEM solver, including an ongoing interactive animation editor. The solver support several elastic materials such as *linear elasticity*, *stvk*, *corotational*, *neohookean* and *fung*. Linear elasticity and stvk are  of special concern to our project. They share same constitutive law but with different strain tensor,
	
	$\Psi(E)=\miu E:E+\frac{\lambda}{2}\text{tr}^2(E)$

thus *small strain tensor* for linear elasticity while *Green strain tensor* for stvk.

## Dynamics ##

I simply use implicit Euler to numerically integrate the Euler-Lagrange equations with Rayleigh damping model, and one can implement other integrator like implicit Newmark as well.
 
## Reduced Solver ##

Up to now, I only implement linear modal analysis in the project. More selection for modal basis such as model direvatives may be incorporated in future. By solving the general eigenvalue problem, we can derive the reduced motion equation. For linear elasticity, the equations reduces to a series of separate second-order ordinary differential equations with constant coefficients, which can be solved analytically given some initial and/or boundary conditions. For stvk materials, we still need to evaluate the reduced force(cubic in reduced coordinates) and the reduced stiffness matrix(quadratic in reduced coordinates), whose coefficients can be pre-computed, see [Barbic05](http://www.ri.cmu.edu/pub_files/pub4/barbic_jernej_2005_1/barbic_jernej_2005_1.pdf) for details.

## Model Warping ##

Linear modal basis will generate articfacts of the simulation results under large deformation. [Huang11](http://www.cad.zju.edu.cn/home/hj/11/dynInterpolation.pdf) proposed Rotation-Strain coordinates to build a bridge connecting the reduced space and the Euclidean space and managed to get the physically plausible results. 

## Animation Editor ##

Ongoing. See [Barbic12](http://www.cs.columbia.edu/cg/pdfs/1344873352-BarbicSinGrinspun-SIGGRAPH2012.pdf).