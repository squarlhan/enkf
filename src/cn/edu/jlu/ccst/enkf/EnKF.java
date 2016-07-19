package cn.edu.jlu.ccst.enkf;

import java.util.Random;

import Jama.Matrix;
/**
 *  Created by Squarlhan  July 2016
	This program uses the ensemble kalman filter to estimate a system's state.
	The state is x_new=f(x,u)+w, where u some input, w the
	Gaussian distributed process noise, and f is a nonlinear function. The measurement 
	is y_new=h(x)+v where h is a nonlinear function and v Gaussian distributed measurement noise.                 


       The algorithm used in this code is referenced from the following:
       S Gillijns et. al., "What Is the Ensemble Kalman Filter and How Well Does it Work?"
       Proceedings of the 2006 American Control Conference,
       Minneapolis, Minnesota, USA, June 14-16, 2006, pp 4448-4453.
 * @author Xiaosong Han
 *
 */
public class EnKF {
	
	private Matrix x_tr;
	private Matrix x_estbar;
	private Matrix ybar;
	
	public Matrix f(Matrix x){
		Matrix re = x;
		re.set(0, 0, x.get(0, 0)+0.1*x.get(1, 0)+0.005);
		re.set(1, 0, x.get(1, 0)+0.1);
		return re;		
	}
	
	public Matrix h(Matrix x){
		Matrix re = new Matrix(1,1);
		re.set(0, 0, x.get(0, 0));
		return re;
	}
	
	public Matrix mean(Matrix x, int flag){
		int r= x.getRowDimension();
		int c = x.getColumnDimension();
		Matrix re;
		if(flag ==2){
			re = new Matrix(r,1);
			for(int i = 0; i<= r-1; i++){
				double sum = 0;
				for(int j = 0; j<=c-1; j++){
					sum+=x.get(i, j);
				}
				re.set(i, 0, sum/c);
			}
		}else{
			re = new Matrix(1,c);
			for(int j = 0; j<=c-1; j++){
				double sum = 0;
				for(int i = 0; i<= r-1; i++){
					sum+=x.get(i, j);
				}
				re.set(0, j, sum/c);
			}
		}
		return re;
	}
	
	public Matrix mergebyr(Matrix a, Matrix b){
		if(a == null){
			return b;
		}
		if(b == null){
			return a;
		}
		int r= a.getRowDimension();
		int c1 = a.getColumnDimension();
		int c2 = b.getColumnDimension();
		Matrix re = new Matrix(r, c1+c2);
		for(int i = 0; i<= r-1; i++){
			for(int j = 0; j<=c1-1; j++){
				re.set(i, j, a.get(i, j));
			}
            for(int j = 0; j<=c2-1; j++){
            	re.set(i, c1+j, b.get(i, j));
			}
		}
		return re;
	}
	
	public void preocess(Matrix x_ini, Matrix z, int num_iterations, Matrix w){
		int dummy = x_ini.getRowDimension();
		int num_members = x_ini.getColumnDimension();
		int p1 = 2;
		int m1 = 1;
		Matrix xvec=null;
		Matrix x_estvec=null;
		Matrix yvec=null;
		Matrix x_est=x_ini.copy();
		Matrix w_n=w.copy();
		Matrix W = new Matrix(p1, num_members);
		Matrix Z = new Matrix(m1, num_members);
		Matrix y = new Matrix(m1, num_members);
		Matrix y_for = new Matrix(m1, num_members);
		Matrix y_forbar = new Matrix(m1, num_members);
		//create measurement noise covariance matrix
		Matrix Zcov = new Matrix(m1, m1);
        for(int j = 0;j<=m1-1;j++){
        	Zcov.set(j, j, Math.pow(z.get(j, j), 2));
		}
		
        for(int i = 0; i<= num_iterations-1; i++){
        	Random rand = new  Random();
        	for(int j = 0;j<=p1-1;j++){
        		w_n.set(j, 0, w_n.get(j, 0)*rand.nextGaussian());
        	}
        	//compute true value of state at next% time step
        	x_tr = f(x_tr).plus(w_n);
        	for(int j = 0;j<=num_members-1;j++){
        		int[] r = {0, dummy-1};
    			int[] c = {j};
        		for(int k = 0;k<=p1-1;k++){
        			//create process noise
            		W.set(k, j, w.get(k, 0)*rand.nextGaussian());
            		//forecast state
            		x_est.set(k, j, f(x_est.getMatrix(r,c)).get(k, 0)+W.get(k, j));
            	}
        		for(int k = 0;k<=m1-1;k++){
        			//create measurement noise
        			Z.set(k, j, z.get(k, 0)*rand.nextGaussian());
        			//make measurement
        			y.set(k, j, h(x_tr).get(k, 0)+Z.get(k, j));
        			//forecast measurement
        			y_for.set(k, j, h(x_est.getMatrix(r,c)).get(k, 0));
        		}       		
        	}
        	
        	x_estbar=mean(x_est,2);
        	ybar=mean(y,2);
        	y_forbar=mean(y_for,2);
        	
        	Matrix Ex = x_est.copy();
        	Matrix Ey = y_for.copy();
        	for(int j = 0;j<=num_members-1;j++){
        		for(int k = 0;k<=p1-1;k++){
        			Ex.set(k, j, x_est.get(k, j)-x_estbar.get(k, 0));
        		}
        		for(int k = 0;k<=m1-1;k++){
        			Ey.set(k, j, y_for.get(k, j)-y_forbar.get(k, 0));
        		}
        	}
        	
        	Matrix Pxy = Ex.times(Ey.transpose()).times(1/(num_members-1));
        	//The addition of Zcov to Pyy is not done in Gillijns et. al but I use it here in case num_members=2 or Pyy is nearly singular
        	Matrix Pyy = Ey.times(Ey.transpose()).times(1/(num_members-1)).plus(Zcov);
        	Matrix K = Pxy.times(Pyy.inverse());
        	x_est.plusEquals(K.times(y.minus(y_for)));
        	
        	xvec = mergebyr(xvec, x_tr);
        	x_estvec=mergebyr(x_estvec, x_estbar);
        	yvec=mergebyr(yvec, ybar);
        }
        x_estbar=mean(x_est,2);
	}
	/**
	 *       Example

      The state and measurement in this example  is taken from Dan Simon, "Kalman Filtering", 
      Embedded Systems Programming,2001.
       

	syms  x1 x2;   %variables must be named x1...xn
	f=[x1+.1*x2+.005;x2+.1];
	h=[x1];
	x_tr=[1;1]; %initial value of state
	x_ini=ones(2,20); %ensemble of initial estimate of the state
	w=[10^-3; .02];  %process noise standard deviation
	z=[10];  %measurement noise standard deviation
	num_iterations=600;
 	num_members=20;
 	[a,b,c]=ensemblekfilter(f,h,x_tr,x_ini,w,z,num_iterations);
	 */
	public void test(){
		x_tr = new Matrix(2,1,1);
		Matrix x_ini = new Matrix(2,20,1);
		int num_iterations=600;
		Matrix w = new Matrix(2,1);
		w.set(0, 0, 0.001);
		w.set(1, 0, 0.02);
		Matrix z = new Matrix(1,1, 10);
		preocess(x_ini, z, num_iterations, w);
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		EnKF en = new EnKF();
		en.test();

	}

}
