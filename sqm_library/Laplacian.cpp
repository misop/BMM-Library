#include "stdafx.h"
#include "Laplacian.h"

void contractMeshGraphCPUCotangent(MeshGraph * pMesh) {
	//filebuf *fb = new filebuf();
	static int logInd = 0;
	logInd++;
	//std::stringstream ss;//create a stringstream
	//ss << "log" << logInd;
	//std::string file = ss.str() + ".txt";
	//fb->open(file, ios::out);
	//ostream *os = new ostream(fb);

	Array2D< double >  L, WL, WH, UPA, LV, O;
	Array2D< double > QR, A, B, V;

	// compute Laplace operator from new points

	float *a = new float[2];
	int *idcs = new int[2];

	L = Array2D< double >(pMesh->numOfVertices, pMesh->numOfVertices, 0.0f);

	for (int i = 0; i < pMesh->numOfVertices; i++)
		for (int j = 0; j < pMesh->numOfVertices; j++)
			if (pMesh->E[i][j] && i != j){

				int ind = 0;
				for (int k = 0; k < pMesh->numOfVertices; k++){
					if (pMesh->E[i][k] && pMesh->E[j][k] && i != k && j != k){
						CVector3 v1 = pMesh->pVerts[k] - pMesh->pVerts[j];
						CVector3 v2 = pMesh->pVerts[k] - pMesh->pVerts[i];

						if (ind > 1) {
							//log(LOG_LEVEL_WARNING, "Vrchol narusajuci 2 manifold: " ,k);
						} else
							a[ind++] = AngleBetweenVectors(v1, v2);
					}
				}
				if (ind == 2){
					//L[i][j] = 1.0f;
					double cota = cotan(a[0]) + cotan(a[1]);
					L[i][j] = cotan(a[0]) + cotan(a[1]);

					if (L[i][j] < -FLT_MAX) {
						L[i][j] = -FLT_MAX;
					} else if (L[i][j] > FLT_MAX) {
						L[i][j] = FLT_MAX;
					}
				} else {
					//log(LOG_LEVEL_WARNING, "ind nieje 2, narusa 2D manifold : ", ind);
				}
			}

			delete[] a;
			a = NULL;

			float maxL = -FLT_MAX;
			float minL = FLT_MAX;

			for (int i = 0; i < pMesh->numOfVertices; i++){
				double sum = 0.0;
				for (int j = 0; j < pMesh->numOfVertices; j++){
					if (pMesh->E[i][j])
						sum -= L[i][j];

					if (pMesh->E[i][j] > maxL)
						maxL = pMesh->E[i][j];

					if (pMesh->E[i][j] < minL)
						minL = pMesh->E[i][j];

				}

				L[i][i] = sum;
			}

			boost::numeric::ublas::matrix<float> Ll(pMesh->numOfVertices,pMesh->numOfVertices);


			for (int i=0; i < pMesh->numOfVertices; i++)
				for (int j=0; j < pMesh->numOfVertices; j++)
					Ll(i,j) = L[i][j];

			//log(0, "Matica L cotangent");
			//log(0, L);

			WL = Array2D< double >(pMesh->numOfVertices, pMesh->numOfVertices, 0.0f);
			WH = Array2D< double >(pMesh->numOfVertices, pMesh->numOfVertices, 0.0f);
			UPA = Array2D< double >(pMesh->numOfVertices, pMesh->numOfVertices, 0.0f);

			for (int i = 0; i < pMesh->numOfVertices; i++) {
				float sum = 0;
				for (int j = 0; j < pMesh->numOfVertices; j++) {
					if (pMesh->E[i][j]) sum += 1;
				}
				WL[i][i] = sum * pMesh->wL;//pMesh->wL;
				WH[i][i] = 3;//pMesh->wH[i];
			}

			//(*os) << "WL" << endl;
			//log(WL, os);
			//(*os) << "WH" << endl;
			//log(WH, os);
			//(*os) << "L" << endl;
			//log(L, os);
			UPA = matmult(WL, L);//WL * L;
			//(*os) << "UPA" << endl;
			//log(UPA, os);

			for (int m = 0; m < UPA.dim1(); m++) {
				for (int n = 0; n < UPA.dim2(); n++)
					if (UPA[m][n] < -FLT_MAX) {
						UPA[m][n] = -FLT_MAX;
						//log(LOG_LEVEL_WARNING,"UPA[m][n] bolo mensie ako -float v sume");
						//log(LOG_LEVEL_WARNING, "UPA[m][n] : " ,UPA[m][n]);
					} else if (UPA[m][n] > FLT_MAX) {
						UPA[m][n] = FLT_MAX;
						//log(LOG_LEVEL_WARNING, "UPA[m][n] bolo vacsie ako float v sume");
						//log(LOG_LEVEL_WARNING, "UPA[m][n] : " ,UPA[m][n]);
					}
			}
			//(*os) << "UPA" << endl;
			//log(UPA, os);
			// now solve the system using TNT and JAMA with QR decomposition

			A = Array2D<double>(2 * pMesh->numOfVertices, pMesh->numOfVertices, 0.0f);

			for (int j = 0; j < pMesh->numOfVertices; j++){
				for (int i = 0; i < pMesh->numOfVertices; i++)
					A[pMesh->numOfVertices + i][j] = (double)WH[i][j];
				for (int i = 0; i < pMesh->numOfVertices; i++)
					A[i][j] = (double)UPA[i][j];
			}
			//(*os) << "A" << endl;
			//log(A, os);

			B = Array2D<double>(2 * pMesh->numOfVertices, 3, 0.0f);

			// X
			for (int i = 0; i < pMesh->numOfVertices; i++){
				CVector3 P = pMesh->pVerts[i];
				float wh = WH[i][i];
				B[pMesh->numOfVertices + i][0] = (double)(pMesh->pVerts[i].x * wh);//pMesh->wH[i]);
				B[pMesh->numOfVertices + i][1] = (double)(pMesh->pVerts[i].y * wh);//pMesh->wH[i]);
				B[pMesh->numOfVertices + i][2] = (double)(pMesh->pVerts[i].z * wh);//pMesh->wH[i]);
			}
			//(*os) << "B" << endl;
			//log(B, os);

			JAMA::QR<double> qr(A);
			V = qr.solve(B);
			//(*os) << "V" << endl;
			//log(V, os);

			//delete[] curOneRingArea;
			//curOneRingArea = NULL;

			for (int i = 0; i < pMesh->numOfVertices; i++)
				pMesh->pVerts[i] = CVector3((float)V[i][0], (float)V[i][1], (float)V[i][2]);


			//fb->close();
			//delete fb;
			//delete os;
}
//---------------------------------------------------------------------------

void computeLaplacian(MeshGraph * pMesh) {
	bool goBack = false;

	float avarageCurOneRingArea = 10.0f;
	float avarageOrigOneRingArea = 0.0f;
	float avarageCurOneRingExtent = 0.0f;
	float avarageOrigOneRingExtent = 0.0f;

	bool isOverThreshold = true;

	pMesh->wL = 0.1;//0.001 * (0.001 * sqrt(avarageCurOneRingArea));

	float lastAvarageCurOneRingExtent = FLT_MAX;

	vector<int*> globalNeighbourhoods;

	//while (*ite <= numOfIter){

	contractMeshGraphCPUCotangent(pMesh);

	pMesh->wL = pMesh->wL * 3.0;
}
/*
void log(TNT::Array2D< double > matrix, ostream *os) {
	//if (level > config.LOG_LEVEL)
	//    return;
	for (int i = 0; i < matrix.dim1(); i++) {
		std::string row;
		for (int j = 0; j < matrix.dim2(); j++) {
			System::String ^ s = System::Convert::ToString((float)matrix[i][j]);
			int floatlength = 5;
			int l = floatlength;
			if (s->Length < floatlength)
				l = s->Length;
			for (int j = 0; j < l; j++)
				row += s[j];
			for (int j = 0; j < floatlength+1-l ; j++)
				row += " ";
		}
		(*os) << row << endl;
	}
	(*os) << endl;
	(*os) << endl;
}

void logB(TNT::Array2D< bool > matrix, ostream *os) {
	//if (level > config.LOG_LEVEL)
	//    return;
	for (int i = 0; i < matrix.dim1(); i++) {
		std::string row;
		for (int j = 0; j < matrix.dim2(); j++) {
			System::String ^ s = System::Convert::ToString((float)matrix[i][j]);
			int floatlength = 4;
			int l = floatlength;
			if (s->Length < floatlength)
				l = s->Length;
			for (int j = 0; j < l; j++)
				row += s[j];
			for (int j = 0; j < floatlength+1-l ; j++)
				row += " ";
		}
		(*os) << row << endl;
	}
	(*os) << endl;
	(*os) << endl;
}
*/