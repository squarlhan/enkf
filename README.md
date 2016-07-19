# enkf
Ensemble Kalman Filter
	This program uses the ensemble kalman filter to estimate a system's state.
	The state is x_new=f(x,u)+w, where u some input, w the
	Gaussian distributed process noise, and f is a nonlinear function. The measurement 
	is y_new=h(x)+v where h is a nonlinear function and v Gaussian distributed measurement noise.                 


       The algorithm used in this code is referenced from the following:
       S Gillijns et. al., "What Is the Ensemble Kalman Filter and How Well Does it Work?"
       Proceedings of the 2006 American Control Conference,
       Minneapolis, Minnesota, USA, June 14-16, 2006, pp 4448-4453.
       
       This Java version is re-writen from a matlab version.
       The m version is https://www.mathworks.com/matlabcentral/fileexchange/31093-ensemble-kalman-filter/content/ensemblekfilter/ensemblekfilter.m.
