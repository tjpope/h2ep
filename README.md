# h2ep

## Hamiltonian 2 Electronic Properties

 The tranmission probability is calculated using the non-equilibrium Green's function method. The hamiltonian from SIESTA is needed in the GOLLUM format. Also, the size of the leads (how many atoms). I assume the leads are symmetric. If this isn't the case for you, go find another code (eg. GOLLUM).                                                       
                                                                     
 The hamiltonian mustr be in this form:     

<img src="https://latex.codecogs.com/svg.latex?\Large&space;H=\begin{pmatrix}XX&XX&XX&XX&XX&XX&XX\\XX&h0&h1&00&00&00&XX\\XX&h1&h0&gl&00&00&XX\\XX&00&gl&hs&gr&00&XX\\XX&00&00&gr&h0&h1&XX\\XX&00&00&00&h1&h0&XX\\XX&XX&XX&XX&XX&XX&XX\end{pmatrix}"title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />              
                                                                     
 You need to input where the first h0 starts and ends.                                            
                                                                     
 This code is parallel and uses the OPENMP enviroment. It also uses LAPACK and BLAS. It also uses stuff from fortran 2003 or later and WILL NOT compile properly with previous versions. It only calculates the transmission probability. If you want the thermopower or the IV characteristics, check out mpft.f - found in the src/ folder.

 Once compiled, for further instructions, type:                         
   
  ./h2ep --help                                               
                                                                     

