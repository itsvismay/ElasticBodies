#include "ImplicitNewmark.h"
#include "globals.h"

// static lbfgsfloatval_t evaluateNewmark(void *impn, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
// 	ImplicitNewmark* in = (ImplicitNewmark*) impn;

// 	double beta =0.25;
// 	//from x to x_k
// 	for(int i=0; i<n; i++){
// 		in->x_k(i) = x[i];
// 	}

// 	in->NewmarkXtoTV(in->x_k, in->TVk);
// 	in->NewmarkCalculateElasticForceGradient(in->TVk, in->forceGradient);
// 	in->NewmarkCalculateForces(in->TVk, in->forceGradient, in->x_k, in->f);

// 	lbfgsfloatval_t fx = 0.0;
// 	for(int i=0; i<n; i++){
// 		fx += 0.5*x[i]*in->massVector(i)*x[i] - in->massVector(i)*in->x_old(i)*x[i] - in->massVector(i)*in->h*in->v_old(i)*x[i] - 0.5*(in->h*in->h)*(1-2*beta)*in->f_old(i)*x[i];
// 		g[i] = in->massVector(i)*x[i] - in->massVector(i)*in->x_old(i) - in->massVector(i)*in->h*in->v_old(i) - (in->h*in->h*0.5)*(1-2*beta)*in->f_old(i) - in->h*in->h*beta*in->f(i);
// 	}

// 	double strainE=0;
// 	for(unsigned int i=0; i<in->M.tets.size(); i++){
// 		strainE += in->M.tets[i].undeformedVol*in->M.tets[i].energyDensity;		
// 	}
// 	fx +=in->h*in->h*beta*strainE;
// 	return fx;
// }

// static int progressNewmark(void *instance,
// 	    const lbfgsfloatval_t *x,
// 	    const lbfgsfloatval_t *g,
// 	    const lbfgsfloatval_t fx,
// 	    const lbfgsfloatval_t xnorm,
// 	    const lbfgsfloatval_t gnorm,
// 	    const lbfgsfloatval_t step,
// 	    int n,
// 	    int k,
// 	    int ls){	

// 	printf("Iteration %d:\n", k);
//     printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
//     printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
//     printf("\n");
//     return 0;	
// }

void ImplicitNewmark::initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT){
	IntegratorAbstract::initializeIntegrator(ph, pM, pTV, pTT);
	ZeroMatrix.resize(3*vertsNum, 3*vertsNum);
	ZeroMatrix.setZero();
	Ident.resize(3*vertsNum, 3*vertsNum);
	Ident.setIdentity();
	forceGradient.resize(3*vertsNum, 3*vertsNum);
	grad_g.resize(3*vertsNum, 3*vertsNum);
	x_k.resize(3*vertsNum);
	v_k.resize(3*vertsNum);
	x_k.setZero();
	v_k.setZero();
	TVk = TV;
}

void ImplicitNewmark::NewmarkTVtoX(VectorXd& x_tv, MatrixXd& TVk){
	x_tv.setZero();
	for(unsigned int i = 0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;

		x_tv(3*indices(0)) = TVk.row(indices(0))[0];
		x_tv(3*indices(0)+1) = TVk.row(indices(0))[1];
		x_tv(3*indices(0)+2) = TVk.row(indices(0))[2];

		x_tv(3*indices(1)) = TVk.row(indices(1))[0];
		x_tv(3*indices(1)+1) = TVk.row(indices(1))[1];
		x_tv(3*indices(1)+2) = TVk.row(indices(1))[2];

		x_tv(3*indices(2)) = TVk.row(indices(2))[0];
		x_tv(3*indices(2)+1) = TVk.row(indices(2))[1];
		x_tv(3*indices(2)+2) = TVk.row(indices(2))[2];

		x_tv(3*indices(3)) = TVk.row(indices(3))[0];
		x_tv(3*indices(3)+1) = TVk.row(indices(3))[1];
		x_tv(3*indices(3)+2) = TVk.row(indices(3))[2];
	}
}

void ImplicitNewmark::NewmarkXtoTV(VectorXd& x_tv, MatrixXd& TVk){
	TVk.setZero();
	for(unsigned int i=0; i < M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		TVk.row(indices(0)) = Vector3d(x_tv(3*indices(0)), x_tv(3*indices(0)+1), x_tv(3*indices(0) +2));
		TVk.row(indices(1)) = Vector3d(x_tv(3*indices(1)), x_tv(3*indices(1)+1), x_tv(3*indices(1) +2));
		TVk.row(indices(2)) = Vector3d(x_tv(3*indices(2)), x_tv(3*indices(2)+1), x_tv(3*indices(2) +2));
		TVk.row(indices(3)) = Vector3d(x_tv(3*indices(3)), x_tv(3*indices(3)+1), x_tv(3*indices(3) +2)); 
	}
	return;
}

void ImplicitNewmark::NewmarkCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient){
	forceGradient.setZero();
	
	vector<Trip> triplets1;
	triplets1.reserve(12*12*M.tets.size());	
	for(unsigned int i=0; i<M.tets.size(); i++){
		//Get P(dxn), dx = [1,0, 0...], then [0,1,0,....], and so on... for all 4 vert's x, y, z
		//P is the compute Force Differentials blackbox fxn

		Vector12d dx(12);
		dx.setZero();
		Vector4i indices = M.tets[i].verticesIndex;
		int kj;
		for(unsigned int j=0; j<12; j++){
			dx(j) = 1;
			MatrixXd dForces = M.tets[i].computeForceDifferentials(TVk, dx);
			kj = j%3;
			//row in order for dfxi/dxi ..dfxi/dzl
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0], dForces(0,0)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+1, dForces(1,0)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[0]+2, dForces(2,0)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1], dForces(0,1)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+1, dForces(1,1)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[1]+2, dForces(2,1)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2], dForces(0,2)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+1, dForces(1,2)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[2]+2, dForces(2,2)));

			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3], dForces(0,3)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+1, dForces(1,3)));
			triplets1.push_back(Trip(3*indices[j/3]+kj, 3*indices[3]+2, dForces(2,3)));
			dx(j) = 0; //ASK check is this efficient?
		}
	}
	forceGradient.setFromTriplets(triplets1.begin(), triplets1.end());
	return;
}

void ImplicitNewmark::NewmarkCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f){
	// //gravity
	f.setZero();
	for(unsigned int i=0; i<M.tets.size(); i++){
		double vertex_mass = M.tets[i].undeformedVol/4;//assume const density 1
		Vector4i indices = M.tets[i].verticesIndex;
		f(3*indices(0)+1) += vertex_mass*gravity;
		f(3*indices(1)+1) += vertex_mass*gravity; 
		f(3*indices(2)+1) += vertex_mass*gravity;
		f(3*indices(3)+1) += vertex_mass*gravity;
	}

	//elastic
	for(unsigned int i=0; i<M.tets.size(); i++){
		Vector4i indices = M.tets[i].verticesIndex;
		MatrixXd F_tet = M.tets[i].computeElasticForces(TVk, simTime%2);
		f.segment<3>(3*indices(0)) += F_tet.col(0);
		f.segment<3>(3*indices(1)) += F_tet.col(1);
		f.segment<3>(3*indices(2)) += F_tet.col(2);
		f.segment<3>(3*indices(3)) += F_tet.col(3);
	}
	// cout<<f<<endl<<endl;
	//damping
	f += rayleighCoeff*forceGradient*(x_k - x_old)/h;
	// cout<<f<<endl<<endl;
	return;
}

void ImplicitNewmark::renderNewtonsMethod(){
	//Implicit Code
	v_k.setZero();
	x_k.setZero();
	x_k = x_old;
	v_k = v_old;
	VectorXd f_old = f;

	forceGradient.setZero();
	bool Nan=false;
	int NEWTON_MAX = 100, i =0;
	double gamma = 0.5;
	double beta =0.25;

	// cout<<"--------"<<simTime<<"-------"<<endl;
	// cout<<"x_k"<<endl;
	// cout<<x_k<<endl<<endl;
	// cout<<"v_k"<<endl;
	// cout<<v_k<<endl<<endl;
	// cout<<"--------------------"<<endl;
	for( i=0; i<NEWTON_MAX; i++){
		grad_g.setZero();
	
		NewmarkXtoTV(x_k, TVk);//TVk value changed in function
		NewmarkCalculateElasticForceGradient(TVk, forceGradient); 
		NewmarkCalculateForces(TVk, forceGradient, x_k, f);

		VectorXd g = x_k - x_old - h*v_old - (h*h/2)*(1-2*beta)*InvMass*f_old - (h*h*beta)*InvMass*f;
		grad_g = Ident - h*h*beta*InvMass*(forceGradient+(rayleighCoeff/h)*forceGradient);

		// VectorXd g = RegMass*x_k - RegMass*x_old - RegMass*h*v_old - (h*h/2)*(1-2*beta)*f_old - (h*h*beta)*f;
		// grad_g = RegMass - h*h*beta*(forceGradient+(rayleighCoeff/h)*forceGradient);
	
		// cout<<"G"<<t<<endl;
		// cout<<g<<endl<<endl;
		// cout<<"G Gradient"<<t<<endl;
		// cout<<grad_g<<endl;

		//solve for delta x
		// Conj Grad
		// ConjugateGradient<SparseMatrix<double>> cg;
		// cg.compute(grad_g);
		// VectorXd deltaX = -1*cg.solve(g);

		// Sparse Cholesky LL^T
		// SimplicialLLT<SparseMatrix<double>> llt;
		// llt.compute(grad_g);
		// VectorXd deltaX = -1* llt.solve(g);

		//Sparse QR 
		SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		sqr.compute(grad_g);
		VectorXd deltaX = -1*sqr.solve(g);

		// CholmodSimplicialLLT<SparseMatrix<double>> cholmodllt;
		// cholmodllt.compute(grad_g);
		// VectorXd deltaX = -cholmodllt.solve(g);
		

		x_k+=deltaX;
		if(x_k != x_k){
			Nan = true;
			break;
		}
		if(g.squaredNorm()<.00000001){
			break;
		}
	}
	if(Nan){
		cout<<"ERROR NEWMARK: Newton's method doesn't converge"<<endl;
		cout<<i<<endl;
		exit(0);
	}
	if(i== NEWTON_MAX){
		cout<<"ERROR NEWMARK: Newton max reached"<<endl;
		cout<<i<<endl;
		exit(0);
	}
	v_old = v_old + h*(1-gamma)*InvMass*f_old + h*gamma*InvMass*f;
	x_old = x_k;
}

// void ImplicitNewmark::renderLBFGS(){
// 	//LBFGS
// 	int N=3*vertsNum;
// 	int i, ret = 0;

// 	double gamma = 0.5;

//     lbfgsfloatval_t fx;
//     lbfgsfloatval_t *x = lbfgs_malloc(N);
//     lbfgs_parameter_t param;
//     if (x == NULL) {
//         printf("ERROR: Newmark Failed to allocate a memory block for variables.\n");
//     }

//     /* Initialize the variables. */
//     x_k.setZero();
//     for (i = 0;i < N; i++) {
       
//        x[i] = x_old(i);
//     }
//     x_k = x_old;
//     v_k = v_old;
//     f_old = f;

//      // Initialize the parameters for the L-BFGS optimization. 
//     lbfgs_parameter_init(&param);

//     /*
//         Start the L-BFGS optimization; this will invoke the callback functions
//         evaluateEuler() and progress() when necessary.
//      */
//     ret = lbfgs(N, x, &fx, evaluateNewmark, progressNewmark, this, &param);
//     if(ret<0){
//     	cout<<"ERROR: Newmark Liblbfgs did not converge, code "<<ret<<endl;
//     	exit(0);
//     }
    
//     v_old = v_old + h*(1-gamma)*InvMass*f_old + h*gamma*InvMass*f;
// 	x_old = x_k;

//     lbfgs_free(x);
// }

void ImplicitNewmark::render(VectorXd& ext_force){
	simTime+=1;
	cout<<"n"<<simTime<<endl;
	IntegratorAbstract::printInfo();
	// renderLBFGS();
	renderNewtonsMethod();

	NewmarkXtoTV(x_old, TV);
	return;
}