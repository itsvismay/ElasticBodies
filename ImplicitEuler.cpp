#include "ImplicitEuler.h"

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> Trip;
typedef Matrix<double, 12, 1> Vector12d;
using namespace alglib;

struct VougaGeometryInfo
{
	Eigen::MatrixX3d restpos;
	Eigen::MatrixX3d oldpos;
	Eigen::MatrixX3d oldvels;
	Eigen::MatrixX4i tets;
	Eigen::VectorXd masses;
	double mu, lambda;
	double h;
	double s;
};

void rep_grad(const real_1d_array &x, double func, void *impe){
	cout<<"*********FUNC***********"<<endl;
	cout<<func<<endl;
	for(int i=0; i< ((ImplicitEuler *)impe)->x_k.rows(); i++)
		cout << x[i] << " ";
	cout << endl;
}

void rep_grad_vouga(const real_1d_array &x, double func, void *impe){
	cout<<"*********FUNC***********"<<endl;
	cout<<func<<endl;
	for(int i=0; i< ((VougaGeometryInfo *)impe)->restpos.rows()*3; i++)
		cout << x[i] << " ";
	cout << endl;
}



void computeEnergyGradient(const VougaGeometryInfo &geom, const Eigen::MatrixX3d &curpos, double &energy, MatrixX3d &gradient)
{
	int ntets = geom.tets.rows();
	int nverts = geom.restpos.rows();
	gradient.resize(nverts, 3);
	gradient.setZero();
	energy = 0;

	for(int i=0; i<ntets; i++)
	{
		Eigen::Matrix3d restdef;
		Eigen::Matrix3d curdef;
		for(int j=0; j<3; j++)
		{
			restdef.col(j) = geom.restpos.row(geom.tets(i,j)).transpose() - geom.restpos.row(geom.tets(i,3)).transpose();
			curdef.col(j) = curpos.row(geom.tets(i,j)).transpose() - curpos.row(geom.tets(i,3)).transpose();
		}
		double restvol = fabs(restdef.determinant())/6.0;

		Eigen::Matrix3d F = curdef*restdef.inverse();
		double J = F.determinant();		
		double I1 = (F.transpose()*F).trace();
		double I1bar = pow(J, -2.0/3.0)*I1;
	
	        energy += restvol*(geom.mu/2.0 * (I1bar - 3.0) + geom.lambda/2.0 * (J-1.0)*(J-1.0));
	
		Eigen::Matrix3d P = geom.mu * (pow(J, -2.0/3.0) * F) - geom.mu/3.0 * I1 * pow(J, -5.0/3.0) * F.inverse().transpose() * J + geom.lambda * (J-1.0) * F.inverse().transpose() * J;
		Eigen::Matrix3d H = restvol*P*restdef.inverse().transpose();
		for(int j=0; j<3; j++)
		{
			gradient.row(geom.tets(i,j)) += H.col(j);
			gradient.row(geom.tets(i,3)) -= H.col(j);
		}
	}
}

void function_grad_vouga(const real_1d_array &y, double &func, real_1d_array &grad, void *impe)
{
	VougaGeometryInfo *geom = (VougaGeometryInfo *)impe;
	int nverts = geom->restpos.rows();
	Eigen::MatrixX3d curpos(nverts, 3);
	for(int i=0; i<nverts; i++)
		for(int j=0; j<3; j++)
		{
			curpos(i,j) = geom->h * y[3*i+j] + geom->h * geom->oldvels(i,j) + geom->oldpos(i,j);			
		}

	Eigen::MatrixX3d grade;
	double energy;
	computeEnergyGradient(*geom, curpos, energy, grade);

	double kineticE = 0.0;
	func = 0;
	for(int i=0; i<nverts; i++)
		for(int j=0; j<3; j++){
			func += 0.5*geom->masses[i]*y[3*i+j]*y[3*i+j]/geom->s;
		}
	kineticE += func;
	// func = 0.0;
	func += energy/geom->s;
	for(int i=0; i<nverts; i++)
		for(int j=0; j<3; j++)
		{
			grad[3*i+j] = geom->masses[i]*y[3*i+j]/geom->s + grade(i,j)*geom->h/geom->s;
		}

	cout<<"-----End of run---------------------"<<endl;
	cout<<"-----"<<"positions"<<endl;
	for(int i=0; i<3*nverts; i++){
		cout<<y[i]<<endl;
	}
	cout<<"kinetic term for func"<<endl;
	cout<<kineticE/geom->s<<endl;
	cout<<"potential term"<<endl;
	cout<<energy/geom->s<<endl;
	// cout<<"kinetic term for g"<<endl;
	// cout<<"Masses"<<endl;
	cout<<"-----"<<func<<endl;
	for(int i=0; i<3*nverts; i++){
		cout<<grad[i]<<endl;
	}
	exit(0);
	
}

int ImplicitEuler::alglibLBFGSVouga(VectorXd& ext_force){
    int N = 3*(vertsNum);
    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
    // using LBFGS method.
    //
    // Several advanced techniques are demonstrated:
    // * upper limit on step size
    // * restart from new point
    //
    real_1d_array x;
    double *positions= new double[N];
   	for(int i=0; i<N; i++){
   		positions[i] = 0;//x_old(i);
   	}

    x.setcontent(N, positions);
    double epsg = sqrt(1e-11)/sqrt(convergence_scaling_paramter);
    cout<<"EPSG"<<endl;
    cout<<epsg<<endl;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0;
    ae_int_t maxits = 0;
    minlbfgsstate state;
    minlbfgsreport rep;
    double teststep = 0;
    // function1_grad(x, teststep, this);
    // cout<<"v-old"<<endl;
    // cout<<v_old<<endl;
    // cout<<"x"<<endl;
    // cout<<x_old<<endl;
    //cout<<"first term"<<endl;
    //cout<<(((x_old+h*v_old) - x_old)/h) - v_old<<endl;
    // x[0]+= h*0.0000001;
    // x[3]+= h*1;
    // x[6]+= h*1;
    // x[9]+= h*1;
    // x[12]+= h*1;
    // first run
    minlbfgscreate(12, x, state);
    minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    //minlbfgssetgradientcheck(state, 1e-8);
    // minlbfgssetstpmax(state, stpmax);
    minlbfgssetxrep(state, true);
    
    VougaGeometryInfo geom;
    geom.restpos.resize(vertsNum,3);
    geom.oldpos.resize(vertsNum, 3);
    geom.oldvels.resize(vertsNum,3);
    geom.masses.resize(vertsNum);
    geom.tets.resize(M.tets.size(), 4);
    for(int i=0; i<vertsNum; i++)
    {
	geom.oldpos.row(i) = x_old.segment<3>(3*i);
	geom.oldvels.row(i) = v_old.segment<3>(3*i);
        geom.masses[i] = RegMass.coeff(3*i,3*i);
    }
    for(int i=0; i<M.tets.size(); i++)
	geom.tets.row(i) =  M.tets[i].verticesIndex;

    // Where is this set? Cheat for now
    geom.restpos = geom.oldpos;
    geom.h = h;
    geom.s = convergence_scaling_paramter;
    geom.mu = M.tets[0].mu;
    geom.lambda = M.tets[0].lambda;
    alglib::minlbfgsoptimize(state, function_grad_vouga, rep_grad_vouga, &geom);
    minlbfgsresults(state, x, rep);

    printf("TERMINATION TYPE: %d\n", int(rep.terminationtype)); // EXPECTED: 4
    cout<<epsg<<endl;
    cout << "final x ";
    for(int i=0; i<N; i++){
    	x_k(i) = x_old[i] + h*v_old[i] + h*x[i];
	cout << x_k(i) <<endl;
    }
    cout << endl;
    v_old = (x_k - x_old)/h;
    x_old = x_k;
    ImplicitXtoTV(x_old, TV);
    return 0;
}

void function1_grad(const real_1d_array &y, double &func, real_1d_array &grad, void *impe) 
{	
	ImplicitEuler* in = (ImplicitEuler*) impe;

    int n = 3*(in->vertsNum - in->fixedVerts.size());

    VectorXd y_vec; y_vec.resize(in->vertsNum*3); y_vec.setZero();
    for(int i=0; i<n; i++){
    	in->x_k(i) = y[i]*in->h + in->h*in->v_old(i)+in->x_old(i);
    	y_vec(i) = y[i];
    }
    // cout<<"Pre bfgs x_k"<<endl;
    // cout<<in->x_k<<endl;

    in->ImplicitXtoTV(in->x_k, in->TVk);//TVk value changed in function
	in->ImplicitCalculateForces(in->TVk, in->forceGradient, in->x_k, in->f);

	//ENERGY SCALAR
	func = 0.0;
	double kineticE =0.0;
	kineticE = (0.5*y_vec.transpose()*in->RegMass*y_vec);
	func += kineticE;
	// func =0.0;
	double strainE=0.0;
	for(unsigned int i=0; i<in->M.tets.size(); i++){
		strainE += in->M.tets[i].energy;		
	}

    func += strainE;

    double gravE =0.0;
    for(unsigned int i=0; i<n/3; i++){
    	gravE += in->massVector(3*i+1)*in->x_k(3*i+1)*gravity;
    }
    func -=gravE;
    func = func/in->convergence_scaling_paramter;
   
    //GRADIENT VECTcout<<"potential term"<<endl;

    VectorXd g = (in->RegMass*y_vec - in->h*in->f)/in->convergence_scaling_paramter;
	for(int i=0; i<n; i++){
		grad[i] = g(i);
	}

	// cout<<"Func"<<endl;
	// cout<<func<<endl;
	// cout<<"GRAD"<<endl;
	// cout<<g<<endl;
	// cout<<"-----End of run---------------------"<<endl;
	// cout<<"-----"<<"Positions"<<endl;
	// cout<<y_vec<<endl;
	// cout<<"kinetic term"<<endl;
	// cout<<kineticE/in->convergence_scaling_paramter<<endl;
	// cout<<"potential term"<<endl;
	// cout<<strainE/in->convergence_scaling_paramter<<endl;
	// cout<<"gravity term"<<endl;
	// cout<<gravE/in->convergence_scaling_paramter<<endl;
	// // cout<<"Masses"<<endl;
	// // cout<<in->RegMass<<endl;
	// cout<<"--**---"<<func<<endl;
	// cout<<"--**---"<<g<<endl;
}

int ImplicitEuler::alglibLBFGSVismay(VectorXd& ext_force){
  	external_f = ext_force;
    int N = 3*(vertsNum - fixedVerts.size());
    // cout<<"N value"<<endl;
    // cout<<N<<endl;
    real_1d_array x;
    double *positions= new double[N];
   	for(int i=0; i<N; i++){
   		positions[i] = 0;//x_old(i);
   	}

   	ImplicitTVtoX(x_k, TV);

    x.setcontent(N, positions);
    double epsg = sqrt(1e-11)/sqrt(convergence_scaling_paramter);
    //  cout<<"EPSG"<<endl;
    // cout<<epsg<<endl;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0;
    ae_int_t maxits = 0;
    minlbfgsstate state;
    minlbfgsreport rep;
    double teststep = 0;

    minlbfgscreate(12, x, state);
    minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    minlbfgssetxrep(state, true);


    alglib::minlbfgsoptimize(state, function1_grad, rep_grad, this);
    minlbfgsresults(state, x, rep);

    // printf("TERMINATION TYPE: %d\n", int(rep.terminationtype)); // EXPECTED: 4
    // cout << "final x ";
    for(int i=0; i<N; i++){
    	x_k(i) = x_old[i] + h*v_old[i] + h*x[i];
		// cout << x_k(i) <<endl;
    }
    cout << endl;
    v_old = (x_k - x_old)/h;
    x_old = x_k;
    ImplicitXtoTV(x_old, TV);
    return 0;
}

void ImplicitEuler::findgBlock(VectorXd& g_block, VectorXd& x, VectorXd& x_old, int ignorePast){
	//VectorXd g = (RegMass*x - RegMass*x_old)/h - RegMass*v_old - h*f;
	VectorXd g = (RegMass*x - RegMass*x_old) - (RegMass*v_old*h) - h*h*f;
	g_block = g.head(ignorePast*3);

}

void ImplicitEuler::renderNewtonsMethod(VectorXd& ext_force){
	//Implicit Code
	v_k.setZero();
	x_k.setZero();
	x_k = x_old;
	x_k(0) += 0.01;//2.1898e-09;
	v_k = v_old;
	clock_t t0 = clock();

	int ignorePastIndex = TV.rows() - fixedVerts.size();
	SparseMatrix<double> forceGradientStaticBlock;
	forceGradientStaticBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);	
	clock_t t1 = clock();
	cout<<"t1 t0"<<endl;
	cout<<"SECONDS"<<double(t1 - t0)/CLOCKS_PER_SEC<<endl;

	SparseMatrix<double> RegMassBlock;
	RegMassBlock.resize(3*ignorePastIndex, 3*ignorePastIndex);
	RegMassBlock = RegMass.block(0, 0, 3*ignorePastIndex, 3*ignorePastIndex);
	clock_t t2 = clock();
	cout<<"t2 t1"<<endl;
	cout<<"SECONDS"<<double(t2 - t1)/CLOCKS_PER_SEC<<endl;

	forceGradient.setZero();
	bool Nan=false;
	int NEWTON_MAX = 100, i =0;
	// cout<<"--------"<<simTime<<"-------"<<endl;
	// cout<<"x_k"<<endl;
	// cout<<x_k<<endl<<endl;
	// cout<<"v_k"<<endl;
	// cout<<v_k<<endl<<endl;
	// cout<<"--------------------"<<endl;
	for( i=0; i<NEWTON_MAX; i++){
		clock_t t3 = clock();
		grad_g.setZero();
		ImplicitXtoTV(x_k, TVk);//TVk value changed in function
		ImplicitCalculateElasticForceGradient(TVk, forceGradient); 
		ImplicitCalculateForces(TVk, forceGradient, x_k, f);
		clock_t t4 = clock();
		cout<<"forces and grad f"<< i<<endl;
		cout<<"SECONDS"<<double(t4 - t3)/CLOCKS_PER_SEC<<endl;
		// VectorXd g_block = x_k - x_old -h*v_old -h*h*InvMass*f;
		// grad_g = Ident - h*h*InvMass*forceGradient - h*rayleighCoeff*InvMass*forceGradient;
		

		//Block forceGrad and f to exclude the fixed verts
		forceGradientStaticBlock = forceGradient.block(0,0, 3*(ignorePastIndex), 3*ignorePastIndex);
		VectorXd g_block;
		findgBlock(g_block, x_k, x_old, ignorePastIndex);
		// VectorXd g = RegMass*x_k - RegMass*x_old - h*RegMass*v_old - h*h*f;
		// VectorXd g_block = g.head(ignorePastIndex*3);
		grad_g = RegMassBlock - h*h*forceGradientStaticBlock - h*rayleighCoeff*forceGradientStaticBlock;

		// Sparse Cholesky LL^T
		if(llt_solver.info() == Eigen::NumericalIssue){
			cout<<"Possibly using a non- pos def matrix in the LLT method"<<endl;
			exit(0);
		}
		clock_t t5 = clock();
		cout<<"grad_g g block"<< i<<endl;
		cout<<"SECONDS - "<<double(t5 - t4)/CLOCKS_PER_SEC<<endl;
		

		// for(unsigned int k = 0; k<grad_g.cols(); k++){
		// 	for(unsigned int j = 0; j<grad_g.rows(); j++){
		// 		cout<<k<<", "<<j<<endl;
		// 		if(grad_g.coeffRef(k, j) != grad_g.coeffRef(k, j)){
		// 			cout<<"______grad g nans__________"<<endl;
		// 			exit(0);
		// 		}
		// 		if(grad_g.coeffRef(k, j) != 0 && CholeskyAnalyze.coeffRef(k, j)==0){
		// 			cout<<"_______Analyze problem_________"<<endl;
		// 			exit(0);
		// 		}
		// 	}
		// }

		llt_solver.factorize(grad_g);
	
		clock_t t6 = clock();
		cout<<"factorize"<< i<<endl;
		cout<<"SECONDS - "<<double(t6 - t5)/CLOCKS_PER_SEC<<endl<<endl;

		
		VectorXd deltaX = -1* llt_solver.solve(g_block);
		
		clock_t t7 = clock();
		cout<<"solver"<< i<<endl;
		cout<<"SECONDS - "<<double(t7 - t5)/CLOCKS_PER_SEC<<endl<<endl;

		x_k.segment(0, 3*(ignorePastIndex)) += deltaX;
		
		//Sparse QR 
		// SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
		// sqr.compute(grad_g);
		// VectorXd deltaX = -1*sqr.solve(g_block);
		// x_k.segment(0, 3*(ignorePastIndex)) += deltaX;

		if(x_k != x_k){
			Nan = true;
			break;
		}

		cout<<"g norm"<<endl;
		cout<<g_block.squaredNorm()/convergence_scaling_paramter<<endl;
		if(g_block.squaredNorm()/convergence_scaling_paramter < sqrt(1e-11)/sqrt(convergence_scaling_paramter)){
			break;
		}	

	}
	if(Nan){
		cout<<"ERROR: Newton's method doesn't converge"<<endl;
		cout<<i<<endl;
		exit(0);
	}
	if(i== NEWTON_MAX){
		cout<<"ERROR: Newton max reached"<<endl;
		cout<<i<<endl;
		exit(0);
	}
	v_old = (x_k - x_old)/h;
	x_old = x_k;
}

void ImplicitEuler::ImplicitCalculateForces( MatrixXd& TVk, SparseMatrix<double>& forceGradient, VectorXd& x_k, VectorXd& f){
	// //gravity
	f.setZero();

	for(unsigned int i=0; i<f.size()/3; i++){
		f(3*i+1) += massVector(3*i+1)*gravity;
	}

	//elastic
	if(solver.compare("newton")==0){
		for(unsigned int i=0; i<M.tets.size(); i++){
			Vector4i indices = M.tets[i].verticesIndex;
			MatrixXd F_tet = M.tets[i].oldComputeElasticForces(TVk, simTime%2);
			f.segment<3>(3*indices(0)) += F_tet.col(0);
			f.segment<3>(3*indices(1)) += F_tet.col(1);
			f.segment<3>(3*indices(2)) += F_tet.col(2);
			f.segment<3>(3*indices(3)) += F_tet.col(3);
		}
	}else{
		for(unsigned int i=0; i<M.tets.size(); i++){
			M.tets[i].computeElasticForces(TVk, f);
		}
	}


	// for(int k=0; k<f.rows(); k++){
	// 	if(fabs(external_f(k))>0.0001){
	// 		f(k) += external_f(k);
	// 	}
	// }
	// f += rayleighCoeff*forceGradient*(x_k - x_old)/h;
	return;
}

void ImplicitEuler::ImplicitCalculateElasticForceGradient(MatrixXd& TVk, SparseMatrix<double>& forceGradient){
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

void ImplicitEuler::render(VectorXd& ext_force){
	simTime+=1;
	cout<<"i"<<simTime<<endl;

	if(solver.compare("newton")==0){
		renderNewtonsMethod(ext_force);

		
	}else if(solver.compare("lbfgsvismay")==0){
		alglibLBFGSVismay(ext_force);
		
	}else if(solver.compare("lbfgsvouga")==0){
		alglibLBFGSVouga(ext_force);
	}else{
		cout<<"Solver not specified properly"<<endl;
		exit(0);
	}

	// cout<<"*******************"<<endl;
	// cout<< "New Pos"<<simTime<<endl;
	// cout<<x_old<<endl<<endl;
	// cout<< "New Vels"<<simTime<<endl;
	// cout<<v_old<<endl;
	// cout<<"*****************"<<endl<<endl;
	IntegratorAbstract::printInfo();

	ImplicitXtoTV(x_old, TV);

}


void ImplicitEuler::initializeIntegrator(double ph, SolidMesh& pM, MatrixXd& pTV, MatrixXi& pTT){
	IntegratorAbstract::initializeIntegrator(ph, pM, pTV, pTT);
	ZeroMatrix.resize(3*vertsNum, 3*vertsNum);
	ZeroMatrix.setZero();
	Ident.resize(3*vertsNum, 3*vertsNum);
	Ident.setIdentity();
	forceGradient.resize(3*vertsNum, 3*vertsNum);
	grad_g.resize(3*vertsNum, 3*vertsNum);
	x_k.resize(3*vertsNum);
	v_k.resize(3*vertsNum);
	external_f.resize(3*vertsNum);
	external_f.setZero();
	x_k.setZero();
	v_k.setZero();
	TVk = TV;
}

void ImplicitEuler::ImplicitXtoTV(VectorXd& x_tv, MatrixXd& TVk){
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

void ImplicitEuler::ImplicitTVtoX(VectorXd& x_tv, MatrixXd& TVk){
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