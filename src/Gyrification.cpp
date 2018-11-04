/*************************************************
*	Gyrification.cpp
*
*	Release: February 2016
*	Update: September 2016
*
*	University of North Carolina at Chapel Hill
*	Department of Computer Science
*	
*	Ilwoo Lyu, ilwoolyu@cs.unc.edu
*************************************************/

#include <algorithm>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include "Gyrification.h"
#include "Util/SurfaceUtil.h"
#include "Util/Geom.h"
#include "Geodesic/Geodesic.h"
#include "Geodesic/GeodesicA.h"

using std::cout;
using std::endl;

Gyrification::Gyrification(void)
{
	m_mesh = NULL;
	m_outer = NULL;
	m_vertex = NULL;
	m_populationArea = 0;
}

Gyrification::~Gyrification(void)
{
	delete m_mesh;
	delete m_outer;
	delete [] m_vertex;
	delete [] m_check_work;
	delete [] m_u;
	delete [] m_pu;
}

void Gyrification::run(double rad)
{
	time_t tstart;
	double elapse;

	cout << "- Geodesic Distance Computation (rad = " << rad << ")" << endl;
	tstart = clock();
	//setupSDist(rad);
	setupGDist(rad);
	elapse = (double)(clock() - tstart) / CLOCKS_PER_SEC;
	cout << elapse << " sec elapsed" << endl;
	
	cout << "- Area Lookup Table" << endl;
	tstart = clock();
	precomputeArea();
	elapse = (double)(clock() - tstart) / CLOCKS_PER_SEC;
	cout << elapse << " sec elapsed" << endl;

	cout << "- Area Adjustment" << endl;
	if (m_populationArea == 0) m_populationArea = totalArea();
	m_adjRatio = totalArea() / m_populationArea;
	m_maxArea *= m_adjRatio;
	m_intv *= m_adjRatio;
	cout << "  Adjusted ratio: " << m_adjRatio << endl;
	cout << "  Adjusted size: " << m_maxArea << endl;
	
	cout << "- Gyrification Index" << endl;
	tstart = clock();
	computeGyrification();
	elapse = (double)(clock() - tstart) / CLOCKS_PER_SEC;
	cout << elapse << " sec elapsed" << endl;
}

void Gyrification::run(const char *map)
{
	time_t tstart;
	double elapse;

	cout << "- Loading: Geodesic Distance" << endl;
	tstart = clock();
	loadGMap(map);
	elapse = (double)(clock() - tstart) / CLOCKS_PER_SEC;
	cout << elapse << " sec elapsed" << endl;
	
	cout << "- Area Lookup Table" << endl;
	tstart = clock();
	precomputeArea();
	elapse = (double)(clock() - tstart) / CLOCKS_PER_SEC;
	cout << elapse << " sec elapsed" << endl;

	cout << "- Area Adjustment" << endl;
	if (m_populationArea == 0) m_populationArea = totalArea();
	m_adjRatio = totalArea() / m_populationArea;
	m_maxArea *= m_adjRatio;
	m_intv *= m_adjRatio;
	cout << "  Adjusted ratio: " << m_adjRatio << endl;
	cout << "  Adjusted size: " << m_maxArea << endl;
	
	cout << "- Gyrification Index" << endl;
	tstart = clock();
	computeGyrification();
	elapse = (double)(clock() - tstart) / CLOCKS_PER_SEC;
	cout << elapse << " sec elapsed" << endl;
}

void Gyrification::setThreads(int nThreads)
{
	m_nThreads = nThreads;
}

void Gyrification::initVertex(void)
{
	int nv = m_mesh->nVertex();
	m_vertex = new point[nv];
	for (int i = 0; i < nv; i++)
	{
		m_vertex[i].type = 0;
		m_vertex[i].group = -1;
		m_vertex[i].id = i;
		m_vertex[i].eid = -1;
		m_vertex[i].adjDist = FLT_MAX;
	}

	m_u = new double*[nv];
	m_pu = new double[nv * 3];
	m_lam1 = new double[nv];
	m_lam2 = new double[nv];
	for (int i = 0; i < nv; i++) m_u[i] = &m_pu[i * 3];
}

void Gyrification::open(const char *mesh, const char *sulcus, const char *gyrus, const char *outer, const char *corr)
{
	// read mesh
	if (m_mesh != NULL) delete m_mesh;
	m_mesh = new Mesh();
	m_mesh->openFile(mesh);
	
	if (m_outer != NULL) delete m_outer;
	m_outer = new Mesh();
	m_outer->openFile(outer);
	
	m_check_work = NULL;
	if (corr != NULL)
	{
		FILE *fp = fopen(corr, "r");
		for (int i = 0; i < m_mesh->nVertex(); i++)
		{
			int id;
			if (fscanf(fp, "%d", &id) == -1)
			{
				cout << "Fatal error: something goes wrong during I/O processing" << endl;
				fclose(fp);
				exit(1);
			}
			m_lookup.push_back(id);
		}
		fclose(fp);
		m_check_work = new bool[m_outer->nVertex()];
	}

	int nv = m_mesh->nVertex();
	initVertex();	
	
	// read sulcus
	if (!strcmp(&sulcus[strlen(sulcus) - 5], ".bary"))
		loadScurveBary(sulcus);
	else
		loadScurve(sulcus);
	
	// read gyrus
	if (!strcmp(&gyrus[strlen(gyrus) - 5], ".bary"))
		loadGcurveBary(gyrus);
	else
		loadGcurve(gyrus);
}

void Gyrification::setKernelInterval(double intv)
{
	m_intv = intv;
}

void Gyrification::setMaxKernelSize(double area)
{
	m_maxArea = area;
}

void Gyrification::setSpeed(double speed)
{
	m_speed = speed;
}

void Gyrification::setupSDist(double distance)
{
	// setup gd
	Geodesic gd(m_mesh);
	const double *dist = gd.dist();
	const int *state = gd.state();

	// closest distances
	vector<int> _start_point;
	for (int i = 0; i < m_mesh->nVertex(); i++)
		if (m_vertex[i].type == 1) _start_point.push_back(i);
	sort(_start_point.begin(), _start_point.end());
	_start_point.erase(unique(_start_point.begin(), _start_point.end()), _start_point.end());

	int *start_point = new int[_start_point.size()];
	for (int i = 0; i < _start_point.size(); i++)
		start_point[i] = _start_point[i];
	gd.perform_front_propagation(start_point, _start_point.size(), NULL, 0);

	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		m_vertex[i].adjDist = dist[i];
	}
}

void Gyrification::setupGDist(double distance)
{
	// setup gd
	//GeodesicA gd(m_mesh);	// only for sulcal xor gyral points
	Geodesic gd(m_mesh);	// for both sulcal xor gyral points
	const double *dist = gd.dist();
	const int *state = gd.state();
	const int *source = gd.source();

	// closest distances
	vector<int> _start_point;
	for (int i = 0; i < m_mesh->nVertex(); i++)
	if (m_vertex[i].type == 1 || m_vertex[i].type == 2) _start_point.push_back(i);	// use both sulcal and gyral points
	sort(_start_point.begin(), _start_point.end());
	_start_point.erase(unique(_start_point.begin(), _start_point.end()), _start_point.end());

	int *start_point = new int[_start_point.size()];
	for (int i = 0; i < _start_point.size(); i++)
		start_point[i] = _start_point[i];
	gd.perform_front_propagation(start_point, _start_point.size(), NULL, 0);

	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		m_vertex[i].adjDist = dist[i];
		m_vertex[i].source = source[i];
	}
}

void Gyrification::setupTensorEqualArea(void)
{
	const double area = 1;
	
	int nv = m_mesh->nVertex();
	T = new double*[nv];
	pT = new double[nv * 9];
	memset(pT, 0, sizeof(double) * nv * 9);
	
	for (int i = 0; i < nv; i++)
	{
		double T0[9];
		
		const int *list = m_mesh->vertex(i)->list();
		int nn = m_mesh->vertex(i)->nNeighbor();
		double lambda[3] = {0.0, FLT_MAX, 0.0};
		
		Vector N = m_mesh->normal(i)->fv();
		Vector V2;

		// gradient
		for (int j = 0; j < nn; j++)
		{
			int id = m_mesh->vertex(i)->list(j);
			Vector E = Vector(m_mesh->vertex(i)->fv(), m_mesh->vertex(id)->fv());
			float p = N * E;
			Vector D = E - N * p;
			double d = D.norm();
			double t = fabs(m_vertex[i].adjDist - m_vertex[id].adjDist) / d;
			
			D.unit();
			
			if (lambda[2] < t)
			{
				lambda[2] = t;
				V2 = D;
			}
		}
		Vector V1 = N.cross(V2); V1.unit();
		for (int j = 0; j < nn; j++)
		{
			int id = m_mesh->vertex(i)->list(j);
			Vector E = Vector(m_mesh->vertex(i)->fv(), m_mesh->vertex(id)->fv());
			float p = N * E;
			Vector D = E - N * p;
			double d = D.norm();
			double t = fabs(m_vertex[i].adjDist - m_vertex[id].adjDist) / d;
			D.unit();
			
			double x = V1 * D * t, y = V2 * D * t;
			if (y == 0)
			{
				lambda[1] = t;
				break;
			}
			if (y == lambda[2]) continue;
			double lambda1 = fabs(x / sqrt(1 - y * y / lambda[2] / lambda[2]));
			if (lambda[1] > lambda1)
				lambda[1] = lambda1;
		}

		if (lambda[1] == 0) lambda[1] = 1e-5;
		double s1 = lambda[1], s2 = lambda[2];
		
		T[i] = &pT[i * 9];
		T[i][0] = s1 * V1[0] * V1[0] + s2 * V2[0] * V2[0];
		T[i][4] = s1 * V1[1] * V1[1] + s2 * V2[1] * V2[1];
		T[i][8] = s1 * V1[2] * V1[2] + s2 * V2[2] * V2[2];
		T[i][1] = T[i][3] = s1 * V1[0] * V1[1] + s2 * V2[0] * V2[1];
		T[i][2] = T[i][6] = s1 * V1[2] * V1[0] + s2 * V2[2] * V2[0];
		T[i][5] = T[i][7] = s1 * V1[1] * V1[2] + s2 * V2[1] * V2[2];
	}
	smoothingTensor(10);
	
	// find the maximum traveltime from the same source
	double *maxDist = new double[nv];
	memset(maxDist, 0, sizeof(double) * nv);
	
	for (int i = 0; i < nv; i++)
	{
		double maxNeighbor = m_vertex[i].adjDist;
		int id = i;
		do
		{
			int maxid = -1;
			int nn = m_mesh->vertex(id)->nNeighbor();
			for (int j = 0; j < nn; j++)
			{
				int id2 = m_mesh->vertex(id)->list(j);
				if (m_vertex[id].source != m_vertex[id2].source) continue;
				if (maxNeighbor < m_vertex[id2].adjDist)
				{
					maxNeighbor = m_vertex[id2].adjDist;
					maxid = id2;
				}
			}
			if (id != -1) maxDist[i] = maxNeighbor;
			id = maxid;
		} while (id != -1);
	}
	
	double *adjParam = new double[nv];
	for (int i = 0; i < nv; i++)
	{
		if (maxDist[i] == 0) maxDist[i] = 1;
		double a = m_speed, b = 1;
		adjParam[i] = (b - a) * (m_vertex[i].adjDist - 0) / (maxDist[i] - 0) + a;
		if (adjParam[i] > 1) adjParam[i] = 1;
	}
	smoothingScalar(adjParam, 10);
	for (int i = 0; i < nv; i++)
		if (adjParam[i] > 1) adjParam[i] = 1;

	for (int i = 0; i < nv; i++)
	{
		// eigenvectors
		double lambda[3];
		double eigv[3][3];
		double T1[3][3] = {{T[i][0], T[i][1], T[i][2]}, {T[i][3], T[i][4], T[i][5]}, {T[i][6], T[i][7], T[i][8]}};
		LinearAlgebra::eig3symmetric(T1, lambda, eigv);
		
		T[i][0] = lambda[0] * eigv[0][0] * eigv[0][0] + lambda[1] * eigv[1][0] * eigv[1][0] + lambda[2] * eigv[2][0] * eigv[2][0];
		T[i][4] = lambda[0] * eigv[0][1] * eigv[0][1] + lambda[1] * eigv[1][1] * eigv[1][1] + lambda[2] * eigv[2][1] * eigv[2][1];
		T[i][8] = lambda[0] * eigv[0][2] * eigv[0][2] + lambda[1] * eigv[1][2] * eigv[1][2] + lambda[2] * eigv[2][2] * eigv[2][2];
		T[i][1] = T[i][3] = lambda[0] * eigv[0][0] * eigv[0][1] + lambda[1] * eigv[1][0] * eigv[1][1] + lambda[2] * eigv[2][0] * eigv[2][1];
		T[i][2] = T[i][6] = lambda[0] * eigv[0][2] * eigv[0][0] + lambda[1] * eigv[1][2] * eigv[1][0] + lambda[2] * eigv[2][2] * eigv[2][0];
		T[i][5] = T[i][7] = lambda[0] * eigv[0][1] * eigv[0][2] + lambda[1] * eigv[1][1] * eigv[1][2] + lambda[2] * eigv[2][1] * eigv[2][2];
		
		// copy tensor info
		memcpy(m_u[i], eigv[1], sizeof(double) * 3);
		m_u[i][0] = eigv[2][0]; m_u[i][1] = eigv[2][1]; m_u[i][2] = eigv[2][2];
		m_lam1[i] = adjParam[i];
		m_lam2[i] = area / m_lam1[i];
	}
	
	delete [] adjParam;
	delete [] maxDist;
}

void Gyrification::computeGyrification(void)
{
	// setup gd
	GeodesicA **gd = new GeodesicA*[m_nThreads];
	const double **dist;
	const int **state;
	const double **area;
	for (int i = 0; i < m_nThreads; i++)
	{
		gd[i] = new GeodesicA(m_mesh);
		dist[i] = gd[i]->dist();
		state[i] = gd[i]->state();
		area[i] = gd[i]->area();
		gd[i]->setArea((const double *)&m_areamap2[0]);	// if outer area is constant
	}
	
	setupTensorEqualArea();
	
	// end_points: gyral points
	cout << "Loading: end points";
	vector<int> _end_gyral_point, _end_sulcal_point;
	for (int i = 0; i < m_mesh->nVertex(); i++)
		if (m_vertex[i].type == 2) _end_gyral_point.push_back(i);
		else if (m_vertex[i].type == 1) _end_sulcal_point.push_back(i);
	sort(_end_gyral_point.begin(), _end_gyral_point.end());
	sort(_end_sulcal_point.begin(), _end_sulcal_point.end());
	_end_gyral_point.erase(unique(_end_gyral_point.begin(), _end_gyral_point.end()), _end_gyral_point.end());
	_end_sulcal_point.erase(unique(_end_sulcal_point.begin(), _end_sulcal_point.end()), _end_sulcal_point.end());
	
	int *end_gyral_point = new int[_end_gyral_point.size()];
	for (int i = 0; i < _end_gyral_point.size(); i++)
		end_gyral_point[i] = _end_gyral_point[i];
	int *end_sulcal_point = new int[_end_sulcal_point.size()];
	for (int i = 0; i < _end_sulcal_point.size(); i++)
		end_sulcal_point[i] = _end_sulcal_point[i];
	cout << endl;
	
	// velocity map
	for (int i = 0; i < m_nThreads; i++)
		gd[i]->setTensor((const double **)m_u, (const double *)m_lam1, (const double *)m_lam2);

	if (m_intv == 0) m_intv = m_maxArea;

	int n = (int)ceil(m_maxArea / m_intv);
	time_t tstart = clock();

	int threadID = 0;
	#pragma omp parallel for private(threadID)
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		if (i % 1000 == 0)
		{
			double elapse = (double)(clock() - tstart) / CLOCKS_PER_SEC;
			cout << "Vertex " << i << ": ";
			cout << elapse << " sec elapsed\n";
			fflush(stdout);
			//tstart = clock();
		}
		gd[threadID]->perform_front_propagation(&i, 1, NULL, 0, 1e9, 0, m_maxArea);

		for (int t = 0; t < n; t++)
		{
			double delta = m_intv * (t + 1);
			if (delta > m_maxArea) delta = m_maxArea;

			m_vertex[i].GI.push_back(kernelArea(&m_vertex[i], dist[threadID], state[threadID], area[threadID], delta));
		}
	}
	cout << endl;

	delete [] end_gyral_point;
	delete [] end_sulcal_point;
	for (int i = 0; i < m_nThreads; i++)
		delete gd[i];
	delete [] gd;
}

void Gyrification::smoothingScalar(double *scalar, int n)
{
	int nv = m_mesh->nVertex();
	double *work = new double[nv];
	for (int iter = 0; iter < n; iter++)
	{
		memset(work, 0, sizeof(double) * nv);
		for (int i = 0; i < nv; i++)
		{
			const int *list = m_mesh->vertex(i)->list();
			int nn = m_mesh->vertex(i)->nNeighbor();
			double wsum = 0;
			
			for (int j = 0; j < nn; j++)
			{
				int id = m_mesh->vertex(i)->list(j);
				Vector E = Vector(m_mesh->vertex(i)->fv(), m_mesh->vertex(id)->fv());
				float w = 1 / E.norm();
				wsum += w;
			
				work[i] += scalar[id] * w;
			}
			work[i] /= wsum;
		}
		memcpy(scalar, work, sizeof(double) * nv);
	}
	delete [] work;
}

void Gyrification::smoothingTensor(int n)
{
	int nv = m_mesh->nVertex();
	double *wTt = new double[nv * 9];
	for (int iter = 0; iter < n; iter++)
	{
		memset(wTt, 0, sizeof(double) * 9 * nv);
		for (int i = 0; i < nv; i++)
		{
			double *wTp = &wTt[i * 9];
			const int *list = m_mesh->vertex(i)->list();
			int nn = m_mesh->vertex(i)->nNeighbor();
			double wsum = 0;
			
			for (int j = 0; j < nn; j++)
			{
				int id = m_mesh->vertex(i)->list(j);
				Vector E = Vector(m_mesh->vertex(i)->fv(), m_mesh->vertex(id)->fv());
				Vector N = m_mesh->normal(i)->fv();
				float p = N * E;
				Vector D = E - N * p;
				float w = 1 / D.norm();
				wsum += w;
			
				for (int k = 0; k < 9; k++)
					wTp[k] += T[id][k] * w;
			}
			for (int j = 0; j < 9; j++)
				wTp[j] /= wsum;
		}
		memcpy(pT, wTt, sizeof(double) * nv * 9);
	}
	delete [] wTt;
}

void Gyrification::precomputeArea(void)
{
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		double area = vertexArea(m_mesh, i);
		m_areamap1.push_back(area);
	}
	int ndefects = 0;
	for (int i = 0; i < m_outer->nVertex(); i++)
	{
		double area = vertexArea(m_outer, i);
		// hull area cannot exceed the original surface
		if (m_areamap1[i] < area)
		{
			ndefects++;
			area = m_areamap1[i];
		}
		m_areamap2.push_back(area);
	}
	cout << "Total " << ndefects << " outer hull regions corrected." << endl;
}

double Gyrification::kernelArea(point *p, const double *dist, const int *state, const double *area, const double size)
{
	double sum1 = 0, sum2 = 0;
	if (m_check_work != NULL)
		memset(m_check_work, 0, sizeof(bool) * m_outer->nVertex());
	double maxdist = 0;
	for (int i = 0; i < m_mesh->nVertex(); i++)
		if (maxdist < dist[i] && dist[i] < 1e6) maxdist = dist[i];

	if (m_check_work != NULL)
	{	
		for (int i = 0; i < m_mesh->nVertex(); i++)
		{
			// for a selected vertex 
			//if (dist[i] < maxdist)
			if (state[i] == GW_GeodesicVertex::kDead)
			{
				if ((area != NULL && size > 0 && area[i] < size) || (area == NULL || size <= 0))
				{
					sum1 += m_areamap1[i];
					if (!m_check_work[m_lookup[i]])
					{
						sum2 += m_areamap2[m_lookup[i]];
						m_check_work[m_lookup[i]] = true;
					}
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < m_mesh->nVertex(); i++)
		{
			// for a selected vertex 
			//if (dist[i] < maxdist)
			if (state[i] == GW_GeodesicVertex::kDead)
			{
				if ((area != NULL && size > 0 && area[i] < size) || (area == NULL || size <= 0))
				{
					sum1 += m_areamap1[i];
					sum2 += m_areamap2[i];
				}
			}
		}
	}
	p->kArea1.push_back(sum1);
	p->kArea2.push_back(sum2);
	return sum1 / sum2;
}

double Gyrification::vertexArea(const Mesh *mesh, int id)
{
	const int *list = mesh->vertex(id)->list();
	int nn = mesh->vertex(id)->nNeighbor();

	double sum = 0;
	for (int j = 0; j < nn; j++)
	{
		int id1 = mesh->vertex(id)->list(j);
		int id2 = mesh->vertex(id)->list((j + 1) % nn);
		Vector E1 = Vector(mesh->vertex(id)->fv(), mesh->vertex(id1)->fv());
		Vector E2 = Vector(mesh->vertex(id)->fv(), mesh->vertex(id2)->fv());
		double a = E1.norm(), b = E2.norm();
		E1.unit(); E2.unit();
		double inner = E1 * E2;
		if (inner > 1) inner = 1;
		if (inner < -1) inner = -1;
		double theta = acos(inner);
		double area = a * b * sin(theta) / 2;
		sum += area / 3;
	}
	
	return sum;
}

void Gyrification::saveGI(const char *filename)
{
	int n = (int)ceil(m_maxArea / m_intv);
	for (int t = 0; t < n; t++)
	{
		char fn[1024];
		double delta = m_intv * (t + 1);
		if (delta > m_maxArea) delta = m_maxArea;
		sprintf(fn, "%s.lgi.map.%d.txt", filename, (int)delta);
		FILE *fp = fopen(fn, "w");
		for (int i = 0; i < m_mesh->nVertex(); i++)
		{
			fprintf(fp, "%f\n", m_vertex[i].GI[t]);
		}
		fclose(fp);
	}
}

void Gyrification::saveKMap1(const char *filename)
{
	int n = (int)ceil(m_maxArea / m_intv);
	for (int t = 0; t < n; t++)
	{
		char fn[1024];
		double delta = m_intv * (t + 1);
		if (delta > m_maxArea) delta = m_maxArea;
		sprintf(fn, "%s.k1.map.%d.txt", filename, (int)delta);
		FILE *fp = fopen(fn, "w");
		for (int i = 0; i < m_mesh->nVertex(); i++)
		{
			fprintf(fp, "%f\n", m_vertex[i].kArea1[t]);
		}
		fclose(fp);
	}
}

void Gyrification::saveKMap2(const char *filename)
{
	int n = (int)ceil(m_maxArea / m_intv);
	for (int t = 0; t < n; t++)
	{
		char fn[1024];
		double delta = m_intv * (t + 1);
		if (delta > m_maxArea) delta = m_maxArea;
		sprintf(fn, "%s.k2.map.%d.txt", filename, (int)delta);
		FILE *fp = fopen(fn, "w");
		for (int i = 0; i < m_mesh->nVertex(); i++)
		{
			fprintf(fp, "%f\n", m_vertex[i].kArea2[t]);
		}
		fclose(fp);
	}
}

void Gyrification::saveGMap(const char *filename)
{
	FILE *fp = fopen(filename, "w");
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		fprintf(fp, "%f %d\n", m_vertex[i].adjDist, m_vertex[i].source);
	}
	fclose(fp);
}

void Gyrification::loadGMap(const char *filename)
{
	int nv = m_mesh->nVertex();
	FILE *fp = fopen(filename, "r");
	for (int i = 0; i < nv; i++)
	{
		float dist;
		int source;
		if (fscanf(fp, "%f %d", &dist, &source) == -1)
		{
			cout << "Fatal error: something goes wrong during I/O processing" << endl;
			fclose(fp);
			exit(1);
		}
		m_vertex[i].adjDist = dist;
		m_vertex[i].source = source;
	}
	fclose(fp);
}

void Gyrification::loadSpoint(const char *filename)
{
	// setup stream
	using std::ifstream;
	char buf[8192];
	ifstream fin;
	int nv = m_mesh->nVertex();

	int group = 0;
	FILE *fp = fopen(filename, "r");
	for (int i = 0; i < nv; i++)
	{
		int flag;
		if (fscanf(fp, "%d", &flag) == -1)
		{
			cout << "Fatal error: something goes wrong during I/O processing" << endl;
			fclose(fp);
			exit(1);
		}
		if (flag == 0) continue;
		int c = 0;
		vector<point> curve;
		{
			m_vertex[i].type = 1;
			m_vertex[i].group = group;
			m_vertex[i].eid = c++;
			curve.push_back(m_vertex[i]);
		}
		if (!curve.empty())
		{
			m_sulcus.push_back(curve);
			group++;
		}
	}
	fclose(fp);
}

void Gyrification::loadScurve(const char *filename)
{
	// setup stream
	using std::ifstream;
	char buf[65536];
	ifstream fin;
	
	int group;
	group = 0;
	fin.open(filename);
	while (fin.getline(buf, sizeof(buf)))
	{
		vector<point> curve;
		char *tokens;
		char *ptr = strtok(buf, " ");
		for (int i = 0; ptr != NULL; i++)
		{
			int id = atoi(ptr);
			if (id != 0 || (id == 0 && *ptr != 13))
			{
//				memcpy(m_vertex[id].tan, umax[id], sizeof(double) * 3);
				m_vertex[id].type = 1;
				m_vertex[id].group = group;
				m_vertex[id].eid = i;
				curve.push_back(m_vertex[id]);
			}
			ptr = strtok(NULL, " \r\n");
		}
		if (!curve.empty())
		{
			m_sulcus.push_back(curve);
			group++;
		}
	}
	fin.close();
}

void Gyrification::loadScurveBary(const char *filename)
{
	// setup stream
	using std::ifstream;
	char buf[65536];
	ifstream fin;
	
	fin.open(filename);
	fin.getline(buf, sizeof(buf));
	int n = atoi(buf);
	for (int group = 0; group < n; group++) 
	{
		vector<point> curve;
		fin.getline(buf, sizeof(buf));
		int p = atoi(buf);
		for (int i = 0; i < p; i++) 
		{
			fin.getline(buf, sizeof(buf));
			int id1, id2;
			float coeff;
			sscanf(buf, "%d %d %f", &id1, &id2, &coeff);
			if (id1 == id2 || coeff == 1)
			{
				m_vertex[id1].type = 1;
				m_vertex[id1].group = group;
				m_vertex[id1].eid = i;
				curve.push_back(m_vertex[id1]);
			}
		}
		if (!curve.empty())
		{
			m_sulcus.push_back(curve);
		}
	}
	fin.close();
}

void Gyrification::loadGcurve(const char *filename)
{
	// setup stream
	using std::ifstream;
	char buf[65536];
	ifstream fin;
	
	int group;
	group = 0;
	fin.open(filename);
	while (fin.getline(buf, sizeof(buf)))
	{
		vector<point> curve;
		char *tokens;
		char *ptr = strtok(buf, " ");
		for (int i = 0; ptr != NULL; i++)
		{
			int id = atoi(ptr);
			if (id != 0 || (id == 0 && *ptr != 13))
			{
				//if (cmax[id] > 0.05)	// curvature threshold
				{
//					memcpy(m_vertex[id].tan, umin[id], sizeof(double) * 3);
					m_vertex[id].type = 2;
					m_vertex[id].group = group;
					m_vertex[id].eid = i;
					curve.push_back(m_vertex[id]);
				}
			}
			ptr = strtok(NULL, " \r\n");
		}
		if (!curve.empty())
		{
			m_gyrus.push_back(curve);
			group++;
		}
	}
	fin.close();
}

void Gyrification::loadGcurveBary(const char *filename)
{
	// setup stream
	using std::ifstream;
	char buf[65536];
	ifstream fin;
	
	fin.open(filename);
	fin.getline(buf, sizeof(buf));
	int n = atoi(buf);
	for (int group = 0; group < n; group++) 
	{
		vector<point> curve;
		fin.getline(buf, sizeof(buf));
		int p = atoi(buf);
		for (int i = 0; i < p; i++) 
		{
			fin.getline(buf, sizeof(buf));
			int id1, id2;
			float coeff;
			sscanf(buf, "%d %d %f", &id1, &id2, &coeff);
			if (id1 == id2 || coeff == 1)
			{
				m_vertex[id1].type = 2;
				m_vertex[id1].group = group;
				m_vertex[id1].eid = i;
				curve.push_back(m_vertex[id1]);
			}
		}
		if (!curve.empty())
		{
			m_sulcus.push_back(curve);
		}
	}
	fin.close();
}

void Gyrification::saveKernelBoundary(const char *filename)
{
	FILE *fp = fopen(filename, "a");
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		for (int j = 0; j < m_vertex[i].boundary.size(); j++)
		{
			fprintf(fp, "%d ", m_vertex[i].boundary[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void Gyrification::saveKernelBoundary(FILE *fp, int index, const double *area, double maxArea)
{
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		bool boundary = false;
		for (int j = 0; j < m_mesh->vertex(i)->nNeighbor() && area[i] <= maxArea && !boundary; j++)
		{
			int id = m_mesh->vertex(i)->list(j);
			if (area[id] > maxArea) boundary = true;
		}
		if (boundary) m_vertex[index].boundary.push_back(i);
	}
	int sz = m_vertex[index].boundary.size();
	fwrite (&sz , sizeof(int), 1, fp);
	for (int i = 0; i < sz; i++)
	{
		int id = m_vertex[index].boundary[i];
		fwrite (&id , sizeof(int), 1, fp);
//		fprintf(fp, "%d ", i);
	}
	m_vertex[index].boundary.clear();
}

void Gyrification::getKernelBoundary(int index, const double *area, double maxArea)
{
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		bool boundary = false;
		for (int j = 0; j < m_mesh->vertex(i)->nNeighbor() && area[i] <= maxArea && !boundary; j++)
		{
			int id = m_mesh->vertex(i)->list(j);
			if (area[id] > maxArea) boundary = true;
		}
		if (boundary) m_vertex[index].boundary.push_back(i);
	}
}

double Gyrification::totalArea(void)
{
	double area = 0;
	for (int i = 0; i < m_mesh->nVertex(); i++)
	{
		area += m_areamap1[i];
	}
	return area;
}

void Gyrification::setPopulationArea(double populationArea)
{
	m_populationArea = populationArea;
}
