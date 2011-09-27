//============================================================================
// Name        : EMTC.cpp
// Author      : Chu, Shumo
// Version     : 0.1
// Description : External Memory Algorithm for Triangle Counting
//============================================================================

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <cstring>
#include <map>
#include <ctime>
#include <vector>
#include "heapbitset.h"
#include "runtimecounter.h"

using namespace std;

//typedef
typedef int vid;
typedef int count_type;	//type of count number

//----------------------const variables---------------------
const long int BLK_SZ=4194304; //4M sizeof blk for I/O
const int SZ_PTR=sizeof(void*);
const int SZ_VID=sizeof(vid);
const int SZ_LONG = sizeof(long);
const int SZ_INT = sizeof(int);
const unsigned vidPerBlk = BLK_SZ/SZ_VID;

//---------------------Experiment Settings ------------------
const int maxNoPart = 164240000; //maximum number of nodes in each partition
const int maxDeg = 4000*1024;	//maxium degree of graph
const int maxNumPart = 300;	//maximum number of partitions |warning this value can not less than 2|
const float scale = 1.5;

//---------------------Global Variables ---------------------
long * tc;
long Gno;
long Gsz;
long maxBlkPart = 32;
count_type * result;  //list of no of local triangles
count_type totalNumber; //total number of triangles
long noRead;
long noWrite;
float partTime;


FILE * iGraph;

char * rbuff;
char * wbuff;

char * mem;
vid ** cGraph;
vid * adjBuff;

//map<vid,int> cmap;
vid * cmap;

//function declaration
void ImCt(int noNodePart);
void output(const char * oName);
void outputGraph(int noNode, const char * gName);

////////////////////////////////////////////////////////
//		sequential partitioning algorithm			  //
////////////////////////////////////////////////////////
void emtcA() {

	FILE * GraphFile = iGraph;
	FILE * tmpGraph;


	int iteration = 0;

	long szCGraph = Gsz;	//size of current graph
	long noCGraph = Gno;	//number of nodes of current graph

	memset(result,0,Gno*sizeof(count_type));

	while (true) {

		tmpGraph = tmpfile();

		printf("begin %d th iteration \n",iteration);
		printf("Graph size: %ld \n",szCGraph);

		long rest = 0;

		int num_part = 0; //number of partition

		long unProcessed = szCGraph;
		long maxNoRead = noCGraph;
		szCGraph = 0;
		noCGraph = 0;

		unsigned ptrWbuff = 0;	//pointer to unused in writing buffer

		long gReadCount = 0;

		//scan the graph once
		while (unProcessed > 0) {

			Runtimecounter rt;
			rt.start();

			int readCount = 0;
			num_part ++;

			printf("processing part: %d \n",num_part);

			int noBlkRead;
			if(unProcessed-rest > maxBlkPart*BLK_SZ ){
				noBlkRead = maxBlkPart;
			}
			else{
				noBlkRead = (unProcessed-rest)/BLK_SZ;
				if((unProcessed-rest)%BLK_SZ!=0) noBlkRead++;
			}

			char * curPos = mem+rest;
			char * endPos = curPos+noBlkRead*BLK_SZ;

			if(endPos>mem+(maxBlkPart+2)*BLK_SZ){
				fprintf(stderr,"block size too small: gReadCount %ld \n",gReadCount);
				exit(1);
			}

			//load one partition to memory
			for(int i=0;i<noBlkRead;i++){
				size_t readNum = fread(curPos,SZ_VID,vidPerBlk,GraphFile);
				if(readNum!=vidPerBlk){
					fprintf(stderr,"readNum %zu \n",readNum);
				}
				noRead++;
				curPos += BLK_SZ;
			}


			vid * adj;	//avoid unnecessary adress redirection

			//cmap.clear();
			memset(cmap,-1,Gno*SZ_VID);

			//initiallize adjacent list
			curPos = mem;

			while(gReadCount<maxNoRead && curPos+2*SZ_VID<endPos){

			   	adj = cGraph[readCount]=(vid *)curPos;
			   	if( (curPos+(2+adj[1])*SZ_VID)>endPos ) break;
				curPos+=(adj[1]+2)*SZ_VID;
				cmap[ adj[0] ] = readCount;
				readCount++;
				gReadCount++;
			}

			rest = endPos-curPos;

			unProcessed -= (curPos-mem);

			if(readCount >= maxNoPart){
				fprintf(stderr,"error, maxNoPart %d readCount %d \n",maxNoPart,readCount);
				exit(1);
			}

			rt.stop();
			partTime+=rt.GetRuntime();

			//count the number of triangles in memory
			//fprintf(stderr,"begin imct \n");
			//ImCt(readCount);
			//fprintf(stderr,"END imct \n");
			//remove inter-partition edges and output part of new graph
			long rmNum = 0;
			for(int i=0;i<readCount;i++){
				adjBuff[0] = cGraph[i][0];
				int i_deg = cGraph[i][1]+2;
				int n_deg = 0;

				for(int j=2;j<i_deg; j++){
					//if(cmap.find( cGraph[i][j] ) == cmap.end()){
					if(cmap[ cGraph[i][j] ] == -1){
						adjBuff[n_deg+2] = cGraph[i][j];
						n_deg ++;

						if(n_deg == maxDeg){
							fprintf(stderr,"adjBuff is too small id %d n_deg %d\n",adjBuff[0],n_deg);
							exit(1);
						}
					}
					else{
						rmNum++;
					}

				}

				if(n_deg>0){

					adjBuff[1] = n_deg;

					//copy contents from adjBuff to writing buffer
					unsigned copyAmt = n_deg+2;

					szCGraph += (n_deg+2)*SZ_VID;
					noCGraph++;

					while(copyAmt > 0){

						if(copyAmt < vidPerBlk-ptrWbuff ){

							memcpy(wbuff+ptrWbuff*SZ_VID,adjBuff+(n_deg+2-copyAmt),copyAmt*SZ_VID);
							ptrWbuff = ptrWbuff+copyAmt;
							copyAmt = 0;
						}
						else{ //copyAmt >= vidPerBlk-ptrWbuff
							memcpy(wbuff+ptrWbuff*SZ_VID,adjBuff+(n_deg+2-copyAmt),(vidPerBlk-ptrWbuff)*SZ_VID);
							size_t writeNum = fwrite(wbuff,SZ_VID,vidPerBlk,tmpGraph);
							if(writeNum != vidPerBlk){
								fprintf(stderr,"writeNum %zu\n",writeNum);
							}

							noWrite++;
							//if(noWrite%10==0)
							//printf("%ld th write \n",noWrite);

							copyAmt -= vidPerBlk-ptrWbuff;
							ptrWbuff = 0;
						}
					}

				}
			}

			memmove(mem,curPos,rest); //move rest to front

			printf("rmNum %ld\n",rmNum);

		}

		//flush writing buffer
		if(ptrWbuff>0){
			fwrite(wbuff,SZ_VID,ptrWbuff,tmpGraph);
			noWrite++;
			ptrWbuff = 0;
		}

		++iteration;

		fclose(GraphFile);

		printf("iteration end, number of partition %d \n",num_part);
		printf("-----------------------------------------\n");

		//exit if the number of partitions is less or equal than 2
		if (num_part <= 2) {
			break;
		}
		else{
			GraphFile = tmpGraph;

			if(GraphFile == NULL){
				fprintf(stderr,"Graph file is null \n");
				exit(1);
			}

			rewind(GraphFile);
		}

		//outputGraph(noCGraph,"newGraph.txt");

	}
	printf("%IT iteration %d\n",iteration);

}


//best effort partitioning
void partB(FILE * bGraph, long szCGraph, int partNum, FILE ** partFile, long * partSz, long * partNo){


	//greedy partitioning
	unsigned ptrRbuff = 0;
	long unProcessed = szCGraph;
	long maxSz = maxBlkPart*BLK_SZ;	//maximum size of a partition

	vid * reBuff= (vid *) rbuff;
	fread(reBuff, SZ_VID, vidPerBlk, bGraph);

	vid * dsPart = (vid *)malloc(Gno*SZ_VID);
	memset(dsPart,-1,Gno*SZ_VID);

	int * counters =  (int *) mem;
	vid * adj = (vid *) mem+partNum*SZ_INT;

	//initialize writing buffers
	char * pbMem = (char *)malloc(partNum*(BLK_SZ+SZ_PTR+SZ_INT));

	vid ** wrBuffs = (vid **)pbMem;
	char * curPos = pbMem+partNum*SZ_PTR;

	for(int i=0;i<partNum;i++){
		wrBuffs[i] = (vid *) curPos;
		curPos += BLK_SZ;
	}

	int * ptrWbuffs = (int *)curPos;
	memset(ptrWbuffs,0,partNum*SZ_INT);

	//partition the graph by scanning the graph once
	int readCount = 0;
	while(unProcessed>0){

		vid u = adj[0] = reBuff[ptrRbuff];

		if((++ptrRbuff)== vidPerBlk){
			fread(reBuff, SZ_VID, vidPerBlk, bGraph);
			ptrRbuff = 0;
		}

		vid u_deg = adj[1] = reBuff[ptrRbuff];

		if ((++ptrRbuff) == vidPerBlk) {
			fread(reBuff, SZ_VID, vidPerBlk, bGraph);
			ptrRbuff = 0;
		}

		//----------decide which partition u belong to------------------
		memset(counters,0,partNum*SZ_VID);	//empty the counters

		for(int i=0; i<u_deg;i++){

			int v = adj[i+2] = reBuff[ptrRbuff];

			if ((++ptrRbuff) == vidPerBlk) {
				fread(reBuff, SZ_VID, vidPerBlk, bGraph);
				ptrRbuff = 0;
			}

			if(dsPart[v]>=0){
				counters[ dsPart[v] ] ++;
			}
		}

		//find the partition with most neigbbors
		int maxCount = -1;
		int maxPart = 0;

		for(int i=0;i<partNum;i++){
			if( counters[i]>maxCount && (partSz[i]+(2+u_deg)*SZ_VID <= maxSz)){
				maxCount = counters[i];
				maxPart = i;
			}
		}

		if(maxCount == -1){
			fprintf(stderr, "rc %d u %d cannot find a valid partition \n",readCount,u);
			exit(1);
		}
		else if(maxCount == 0){
			int minSz = partSz[0];
			for (int i = 0; i < partNum; i++) {
				if (partSz[i] < minSz && (partSz[i] + (2 + u_deg) * SZ_VID <= maxSz)) {
					maxPart = i;
					minSz = partSz[i];
				}
			}
		}

		//copy u and its adjacent list to the buffer
		int cpAmt = u_deg+2;
		unProcessed -= cpAmt*SZ_VID;
		partSz[maxPart] += cpAmt*SZ_VID;
		partNo[maxPart] ++;
		dsPart[u] = maxPart;

		while(cpAmt>0){
			if ((unsigned)cpAmt < vidPerBlk-ptrWbuffs[maxPart]) {

				memcpy(wrBuffs[maxPart]+ptrWbuffs[maxPart],adj+(u_deg+2-cpAmt),cpAmt*SZ_VID);
				ptrWbuffs[maxPart] += cpAmt;
				cpAmt = 0;

			} else {//cpAmt > vidPerBlk-ptrWbuffs[maxPart]

				memcpy(wrBuffs[maxPart]+ptrWbuffs[maxPart],adj+(u_deg+2-cpAmt),(vidPerBlk-ptrWbuffs[maxPart])*SZ_VID);
				fwrite(wrBuffs[maxPart],SZ_VID,vidPerBlk,partFile[maxPart]);
				cpAmt -= vidPerBlk-ptrWbuffs[maxPart];
				ptrWbuffs[maxPart] = 0;
				noWrite++;
			}
		}

		readCount++;
	}

	//flush buffer
	for(int i=0;i<partNum;i++){
		if(ptrWbuffs[i]>0){
			fwrite(wrBuffs[i],SZ_VID,ptrWbuffs[i],partFile[i]);
			noWrite++;
		}
	}

	//output information
	delete dsPart;
	delete pbMem;

}

//triangle counting algorithm using best effort partitioning
void emtcB(){

	long szCGraph = Gsz;	//initial size of current graph
	long noCGraph = Gno;	//initial number of nodes of current graph

	FILE * GraphFile = iGraph;
	FILE * tmpGraph;	//temp file
	FILE ** partFile = (FILE **)malloc(maxNumPart*SZ_PTR); //partition files

	long * partSz = (long *)malloc(maxNumPart*SZ_LONG);
	long * partNo = (long *)malloc(maxNumPart*SZ_LONG);

	int iteration = 0;

	while(true){

		printf("begin %d th iteration \n",iteration);
		printf("Graph size: %ld \n",szCGraph);

		tmpGraph = tmpfile();

		//decide number of partitions
		int num_part;

		if( szCGraph > (maxBlkPart*BLK_SZ)*2 ){
			num_part = szCGraph/(maxBlkPart*BLK_SZ)+1;
		}
		else{
			num_part = 2;
		}

		//initialize files
		for(int i=0;i<num_part;i++){
			partFile[i] = tmpfile();
			partSz[i] = 0;
			partNo[i] = 0;
		}

		rewind(GraphFile);
		partB(GraphFile, szCGraph, num_part, partFile,partSz,partNo);

		long ptrWbuff = 0;

		szCGraph = 0;
		noCGraph = 0;

		for(int i = 0;i<num_part; i++){

			printf("processing part: %d \n",i);

			long cPartNo = partNo[i];	//get number of vertices of part i

			int numBlkRead,readCount = 0;

			partSz[i]%BLK_SZ==0?numBlkRead = partSz[i]/BLK_SZ:numBlkRead = partSz[i]/BLK_SZ+1;

			char * curPos = mem;
			char * endPos = mem+numBlkRead*BLK_SZ;

			vid * adj;	//avoid unnecessary adress redirection

			//load partition i into memory

			rewind(partFile[i]);
			//cmap.clear();
			memset(cmap,-1,Gno*SZ_VID);

			for(int j=0;j<numBlkRead;j++){
				size_t readNum = fread(curPos, SZ_VID, vidPerBlk, partFile[i]);
				if (readNum != (size_t) vidPerBlk) {
					fprintf(stderr, "readNum %zu \n", readNum);
				}

				noRead++;
				curPos += BLK_SZ;
			}

			curPos = mem;

			while (readCount < cPartNo) {

				adj = cGraph[readCount] = (vid *) curPos;

				if ((curPos + (2 + adj[1]) * SZ_VID) > endPos){
					fprintf(stderr,"error,rc %d part[%d] is out of bound\n",readCount,i);
					exit(1);
				}

				curPos += (adj[1] + 2) * SZ_VID;
				cmap[ adj[0] ] = readCount;
				readCount++;

			}

			if(readCount!=cPartNo){
				fprintf(stderr,"readCount %d cPartNo %ld",readCount,cPartNo);
				exit(1);
			}

			//in memory triangle computation
			ImCt(readCount);

			//remove intra-partition edges and write new graph file
			for (int i = 0; i < readCount; i++) {
				adjBuff[0] = cGraph[i][0];
				vid i_deg = cGraph[i][1] + 2;
				vid n_deg = 0;

				for (int j = 2; j < i_deg; j++) {
					//if (cmap.find(cGraph[i][j]) == cmap.end()) {
					if (cmap[ cGraph[i][j] ] == -1) {

						adjBuff[n_deg + 2] = cGraph[i][j];
						n_deg++;

						if (n_deg == maxDeg) {
							fprintf(stderr,"adjBuff is too small id %d n_deg %d\n",adjBuff[0], n_deg);
							exit(1);
						}
					}
				}

				if (n_deg > 0) {

					adjBuff[1] = n_deg;

					//copy contents from adjBuff to writing buffer
					int copyAmt = n_deg+2;

					szCGraph += (n_deg + 2)*SZ_VID;
					noCGraph++;

					while (copyAmt > 0) {

						if (copyAmt < vidPerBlk - ptrWbuff) {

							memcpy(wbuff + ptrWbuff * SZ_VID, adjBuff + (n_deg
									+ 2 - copyAmt), copyAmt * SZ_VID);
							ptrWbuff = ptrWbuff + copyAmt;
							copyAmt = 0;

						} else { //copyAmt >= vidPerBlk-ptrWbuff

							memcpy(wbuff + ptrWbuff * SZ_VID, adjBuff + (n_deg
									+ 2 - copyAmt), (vidPerBlk - ptrWbuff)
									* SZ_VID);
							size_t writeNum = fwrite(wbuff, SZ_VID, vidPerBlk,
									tmpGraph);
							if (writeNum != (size_t)vidPerBlk) {
								fprintf(stderr, "writeNum %zu\n", writeNum);
							}

							noWrite++;
							//if (noWrite % 10 == 0)
							//	printf("%ld th write \n", noWrite);

							copyAmt -= vidPerBlk - ptrWbuff;
							ptrWbuff = 0;

						}
					}

				}
			}


		}

		//flush buffer
		if(ptrWbuff>0){
			fwrite(wbuff,SZ_VID,ptrWbuff,tmpGraph);
			ptrWbuff = 0;
		}

		for(int i=0;i<num_part;i++){
			fclose( partFile[i] );
		}

		printf("iteration end, number of partition %d \n",num_part);
		printf("-----------------------------------------\n");

		if(num_part == 2){
			break;
		}
		else{
			fclose(GraphFile);
			GraphFile = tmpGraph;
		}

		iteration++;
	}

	output("LJtri.txt");

	delete partFile;
	delete partSz;
	delete partNo;


}



//dominating set based partitioning
void partC(FILE * cGraph, long noCGraph, long szCGraph, int partNum, FILE ** partFile, long * partSz, long * partNo){


	Runtimecounter rt;
	rt.start();


	//printf("noCGraph %ld\n",noCGraph);
	//-------------computing dominating set-----------------
	HeapBitset bitmap(Gno);
	bitmap.reset();

	unsigned ptrRbuff = 0;	//pointer for read buffer
	long maxSz = maxBlkPart*BLK_SZ;	//max size of each partition

	short * dsPart = (short *)malloc(Gno*sizeof(short)); //list of partition number of each vertex in DS
	//int * dsList = (int *)malloc(Gno*SZ_INT); //list of vertex in dominating set
	vector<int> dsList;
	int dsCount = 0;

	memset(dsPart,-1,Gno*sizeof(short));
	//memset(dsList,0,Gno*SZ_INT);

	int * counters =  (int *) mem;
	vid * adj = (vid *) mem+partNum*SZ_INT;

	vid * reBuff= (vid *) rbuff;
	fread(reBuff, SZ_VID, vidPerBlk, cGraph);

	//initialize writing buffers
	char * pbMem = (char *) malloc(partNum * (BLK_SZ + SZ_PTR + SZ_INT));

	vid ** wrBuffs = (vid **) pbMem;
	char * curPos = pbMem + partNum * SZ_PTR;

	for (int i = 0; i < partNum; i++) {
		wrBuffs[i] = (vid *) curPos;
		curPos += BLK_SZ;
	}

	int * ptrWbuffs = (int *) curPos;
	memset(ptrWbuffs, 0, partNum * SZ_INT);

	for(long i=0;i<noCGraph;i++){

		vid u = adj[0] = reBuff[ptrRbuff];

		//checking validity of u here
		if(u>Gno || u<0){
			fprintf(stderr,"partc u: %d i %d\n",u,i);
			exit(1);
		}

		if((++ptrRbuff)== vidPerBlk){
			fread(reBuff, SZ_VID, vidPerBlk, cGraph);
			ptrRbuff = 0;
		}

		vid u_deg = adj[1] = reBuff[ptrRbuff];

		if ((++ptrRbuff) == vidPerBlk){
			fread(reBuff, SZ_VID, vidPerBlk, cGraph);
			ptrRbuff = 0;
		}

		if (bitmap[u] == false) {

			bitmap[u] = true;

			for (int j = 0; j < u_deg; j++) {

				vid v = adj[j + 2] = reBuff[ptrRbuff];

				if ((++ptrRbuff) == vidPerBlk) {
					fread(reBuff, SZ_VID, vidPerBlk, cGraph);
					ptrRbuff = 0;
				}

				bitmap[v] = true;
			}

			//add u into DS
			//dsList[dsCount] = u;
			dsList.push_back(u);
			dsCount++;

			//find a partition of u
			int cPart = rand()%partNum;
			if(partSz[cPart]+u_deg+2>maxSz){
				bool flag = true;
				for(int j=1;j<partNum;j++){
					if(partSz[ (cPart+j)%partNum ]+u_deg+2<=maxSz){
						cPart = (cPart+j)%partNum;
						flag = false;
						break;
					}
				}
				if(flag){
					fprintf(stderr,"error, can not find a partition \n");
					exit(1);
				}
			}

			//add u into this partition
			dsPart[u] = (short) cPart;

			int cpAmt = u_deg+2;
			partSz[cPart] += cpAmt*SZ_VID;
			partNo[cPart] ++;

			while(cpAmt>0){
				if ((unsigned)cpAmt < vidPerBlk-ptrWbuffs[cPart]) {

					memcpy(wrBuffs[cPart]+ptrWbuffs[cPart],adj+(u_deg+2-cpAmt),cpAmt*SZ_VID);
					ptrWbuffs[cPart] += cpAmt;
					cpAmt = 0;

				} else {//cpAmt > vidPerBlk-ptrWbuffs[maxPart]

					memcpy(wrBuffs[cPart]+ptrWbuffs[cPart],adj+(u_deg+2-cpAmt),(vidPerBlk-ptrWbuffs[cPart])*SZ_VID);
					fwrite(wrBuffs[cPart],SZ_VID,vidPerBlk,partFile[cPart]);
					cpAmt -= vidPerBlk-ptrWbuffs[cPart];
					ptrWbuffs[cPart] = 0;
				}
			}

		}
		else{//if bitmap[u] != false
			for (int j = 0; j < u_deg; j++) {
				if ((++ptrRbuff) == vidPerBlk) {
					fread(reBuff, SZ_VID, vidPerBlk, cGraph);
					ptrRbuff = 0;
				}
			}
		}
	}

	printf("DS size %d\n",dsCount);

	//-------------partitioning based on DS------------------

	rewind(cGraph);
	fread(reBuff, SZ_VID, vidPerBlk, cGraph);
	ptrRbuff = 0;

	for(long i=0;i<noCGraph;i++){

		vid u = adj[0] = reBuff[ptrRbuff];

		if ((++ptrRbuff) == vidPerBlk) {
			fread(reBuff, SZ_VID, vidPerBlk, cGraph);
			ptrRbuff = 0;
		}

		vid u_deg = adj[1] = reBuff[ptrRbuff];

		if ((++ptrRbuff) == vidPerBlk) {
			fread(reBuff, SZ_VID, vidPerBlk, cGraph);
			ptrRbuff = 0;
		}

		int cpAmt = u_deg+2;

		if(dsPart[u]>=0){//if u is in DS
			cpAmt = u_deg;
			while(cpAmt>0){
				if(ptrRbuff+cpAmt<vidPerBlk){
					ptrRbuff+=cpAmt;
					cpAmt = 0;
				}
				else{
					cpAmt -= (vidPerBlk-ptrRbuff);
					fread(reBuff, SZ_VID, vidPerBlk, cGraph);
					ptrRbuff = 0;
				}
			}
		}
		else{//if u is not in DS

			memset(counters,0,partNum*SZ_INT);

			for (int j = 0; j < u_deg; j++) {

				vid v = adj[j + 2] = reBuff[ptrRbuff];

				if ((++ptrRbuff) == vidPerBlk) {
					fread(reBuff, SZ_VID, vidPerBlk, cGraph);
					ptrRbuff = 0;
				}

				if(dsPart[v]>=0){
					counters[dsPart[v] ]++;
				}
			}

			//decide which partition u should be added
			int maxCount = -1;
			int maxPart = 0;

			for(int j=0;j<partNum;j++){
				if( counters[j]>maxCount && (partSz[j]+(2+u_deg)*SZ_VID <= maxSz)){
					maxCount = counters[j];
					maxPart = j;
				}
			}

			if(maxCount == -1){
				for(int j=0;j<partNum;j++){
					if( counters[j]>maxCount && (partSz[j]+(2+u_deg)*SZ_VID <= maxSz+BLK_SZ)){
						maxCount = counters[j];
						maxPart = j;
					}
				}

				if(maxCount == -1){
					fprintf(stderr, "cannot find a valid partition u_deg %d\n",u_deg);
					for(int j = 0;j < partNum;j++){
						fprintf(stderr,"j: %d  partsz: %ld\n",j,partSz[j]);
					}
					exit(1);
				}
			}

			//copy u and its adjacent list to the buffer
			int cpAmt = u_deg+2;
			partSz[maxPart] += cpAmt*SZ_VID;
			partNo[maxPart] ++;
			dsPart[u] = maxPart;

			while(cpAmt>0){
				if ((unsigned)cpAmt < vidPerBlk-ptrWbuffs[maxPart]) {

					memcpy(wrBuffs[maxPart]+ptrWbuffs[maxPart],adj+(u_deg+2-cpAmt),cpAmt*SZ_VID);
					ptrWbuffs[maxPart] += cpAmt;
					cpAmt = 0;

				} else {//cpAmt > vidPerBlk-ptrWbuffs[maxPart]

					memcpy(wrBuffs[maxPart]+ptrWbuffs[maxPart],adj+(u_deg+2-cpAmt),(vidPerBlk-ptrWbuffs[maxPart])*SZ_VID);
					fwrite(wrBuffs[maxPart],SZ_VID,vidPerBlk,partFile[maxPart]);
					cpAmt -= vidPerBlk-ptrWbuffs[maxPart];
					ptrWbuffs[maxPart] = 0;
				}
			}
		}
	}

	//flush buffers
	for(int i=0;i<partNum;i++){
		if(ptrWbuffs[i]>0){
			fwrite(wrBuffs[i],SZ_VID,ptrWbuffs[i],partFile[i]);
		}
	}

	delete pbMem;
	delete dsPart;

	rt.stop();
	partTime+=rt.GetRuntime();

	//delete dsList;
}

//triangle counting algorithm using dominating set based partitioning
void emtcC(){

	long szCGraph = Gsz; //initial size of current graph
	long noCGraph = Gno; //initial number of nodes of current graph

	FILE * GraphFile = iGraph;
	FILE * tmpGraph; //temp file
	FILE ** partFile = (FILE **) malloc(maxNumPart * SZ_PTR); //partition files

	long * partSz = (long *) malloc(SZ_LONG * maxNumPart);	//list of partition size
	long * partNo = (long *) malloc(SZ_LONG * maxNumPart);	//list of node No. of each partition

	int iteration = 0;

	while(true){

		printf("begin %d th iteration \n",iteration);
		printf("Graph size: %ld \n",szCGraph);

		//------------determine number of partitions---------------------
		int num_part;

		if( szCGraph > (maxBlkPart*BLK_SZ)*2 ){
			num_part = szCGraph/(maxBlkPart*BLK_SZ)+1;
		}
		else if(szCGraph < maxBlkPart*BLK_SZ){
			num_part = 1;
		}
		else{
			num_part = 2;
		}

		memset(partSz,0,SZ_LONG * maxNumPart);
		memset(partNo,0,SZ_LONG * maxNumPart);

		tmpGraph = tmpfile();

		//call partitioning algorithm
		for(int i=0;i<num_part;i++){
			partFile[i] = tmpfile();
			partSz[i] = 0;
			partNo[i] = 0;
		}

		rewind(GraphFile);
		if(num_part > 1){
			partC(GraphFile, noCGraph, szCGraph, num_part, partFile,partSz,partNo);
		}
		else{
			partFile[0] = GraphFile;
			partNo[0] = noCGraph;
			partSz[0] = szCGraph;
		}
		long ptrWbuff = 0;

		szCGraph = 0;
		noCGraph = 0;

		for(int i = 0;i<num_part; i++){

			printf("processing part: %d \n",i);

			long cPartNo = partNo[i];	//get number of vertices of part i

			int numBlkRead,readCount = 0;

			partSz[i]%BLK_SZ==0?numBlkRead = partSz[i]/BLK_SZ:numBlkRead = partSz[i]/BLK_SZ+1;

			char * curPos = mem;
			char * endPos = mem+numBlkRead*BLK_SZ;

			vid * adj;	//avoid unnecessary adress redirection

			//load partition i into memory

			rewind(partFile[i]);
			//cmap.clear();
			memset(cmap,-1,Gno*SZ_VID);

			for(int j=0;j<numBlkRead;j++){
				size_t readNum = fread(curPos, SZ_VID, vidPerBlk, partFile[i]);
				if (readNum != (size_t) vidPerBlk) {
					fprintf(stderr, "readNum %zu \n", readNum);
				}

				noRead++;
				curPos += BLK_SZ;
			}

			curPos = mem;

			while (readCount < cPartNo) {

				adj = cGraph[readCount] = (vid *) curPos;

				if ((curPos + (2 + adj[1]) * SZ_VID) > endPos){
					fprintf(stderr,"2error, rc %d part[%d] out bound, nbr %d deg %d p \n",readCount,i,numBlkRead,adj[0]);
					exit(1);
				}

				curPos += (adj[1] + 2) * SZ_VID;
				cmap[ adj[0] ] = readCount;
				readCount++;

			}

			if(readCount!=cPartNo){
				fprintf(stderr,"readCount %d cPartNo %ld",readCount,cPartNo);
				exit(1);
			}


			//in memory triangle computation
			//ImCt(partNo[i]);

			//remove intra-partition edges and write new graph file
			for (int k = 0; k < readCount; k++) {
				adjBuff[0] = cGraph[k][0];
				vid i_deg = cGraph[k][1] + 2;
				vid n_deg = 0;

				for (int j = 2; j < i_deg; j++) {

					if (cmap[ cGraph[k][j] ] == -1) {

						adjBuff[n_deg + 2] = cGraph[k][j];
						n_deg++;

						if (n_deg == maxDeg) {
							fprintf(stderr,"adjBuff is too small id %d n_deg %d\n",adjBuff[0], n_deg);
							exit(1);
						}
					}
				}

				if (n_deg > 0) {

					adjBuff[1] = n_deg;

					//copy contents from adjBuff to writing buffer
					int copyAmt = n_deg+2;

					szCGraph += (n_deg + 2)*SZ_VID;
					noCGraph++;

					while (copyAmt > 0) {

						if (copyAmt < vidPerBlk - ptrWbuff) {

							memcpy(wbuff + ptrWbuff * SZ_VID, adjBuff + (n_deg
									+ 2 - copyAmt), copyAmt * SZ_VID);
							ptrWbuff = ptrWbuff + copyAmt;
							copyAmt = 0;

						} else { //copyAmt >= vidPerBlk-ptrWbuff

							memcpy(wbuff + ptrWbuff * SZ_VID, adjBuff + (n_deg
									+ 2 - copyAmt), (vidPerBlk - ptrWbuff)
									* SZ_VID);
							size_t writeNum = fwrite(wbuff, SZ_VID, vidPerBlk,
									tmpGraph);
							if (writeNum != (size_t)vidPerBlk) {
								fprintf(stderr, "writeNum %zu\n", writeNum);
							}

							noWrite++;
							//if (noWrite % 10 == 0)
							//	printf("%ld th write \n", noWrite);

							copyAmt -= vidPerBlk - ptrWbuff;
							ptrWbuff = 0;

						}
					}

				}
			}


		}//end for

		//flush buffer
		if(ptrWbuff>0){
			fwrite(wbuff,SZ_VID,ptrWbuff,tmpGraph);
			ptrWbuff = 0;
		}

		for(int i=0;i<num_part;i++){
			fclose(partFile[i]);
		}

		printf("-----------------------------------------------\n");


		if(num_part <= 2){
			break;
		}
		else{
			fclose(GraphFile);
			GraphFile = tmpGraph;
		}

		iteration++;
	}//end while

	//output the result
	//output("LJtri.txt");
	printf("%IT iteration %d\n",iteration);

	free(partFile);
	free(partSz);
	free(partNo);
}

//count the number of triangles in memory
void ImCt(int noNodePart){


	fprintf(stderr,"begin no.node %d\n",noNodePart);
	HeapBitset bitmap(Gno);
	bitmap.reset();

	for(int i=0;i<noNodePart;i++){

		if(i%10000 == 0){
			fprintf(stderr,"\r i %d",i);
		}
		int v = cGraph[i][0];
		int v_deg = cGraph[i][1]+2;
		count_type tv = 0;

		int j = 2;

		int startJ = 0;

		while(j < v_deg){
			bitmap[ cGraph[i][j] ] = true;
			if( startJ==0 && cGraph[i][j] > v){//look for start point
				startJ = j;
			}
			j++;
		}

		if(startJ!=0){
			j = startJ;
		}

		while( j < v_deg){

			vid u = cGraph[i][j];
			count_type tu = 0;
			vid it = cmap[u];

			if(it>=0){
				//int idx = (*it).second;
				vid idx = it;
				int u_deg = cGraph[idx][1]+2;


				for(int k=2;k<u_deg;k++){

					if (bitmap[ cGraph[idx][k] ] == true) {

						if ( cGraph[idx][k] > u) {
							result[ cGraph[idx][k] ]++;
							tu++;
							tv++;
						} else {
							//map<vid, int>::iterator iti = cmap.find(cGraph[idx][k]);
							//if (iti == cmap.end()) {
							if(cmap[ cGraph[idx][k] ]==-1){
								result[ cGraph[idx][k] ]++;
								tu++;
								tv++;
							}
						}
					}
				}
			}

			result[u] += tu;
			j++;
		}

		result[v] += tv;

		//recover bitmap
		for(int k=2;k<v_deg;k++){
			bitmap[ cGraph[i][k] ] = false;
		}


	}

}

//output counting result
void output(const char * oName){

	FILE * outfile = fopen(oName,"w");

	totalNumber = 0;

	for( int i = 0; i < Gno; i++){
		fprintf(outfile,"%d:%ld\n",i,result[i]);
	}

	fclose(outfile);
}

//output graph in memory
void outputGraph(int noNode, const char * gName){

	FILE * ogFile = fopen(gName,"w");

	if(noNode > Gno && noNode<0){
		fprintf(stderr,"noNode is wrong %d\n",noNode);
		exit(1);
	}

	if(ogFile == NULL){
		fprintf(stderr,"cannot open file %s \n",gName);
		exit(1);
	}



	for(int i = 0; i < noNode; i++){

		int deg = cGraph[i][1];
		fprintf(ogFile,"%d,%d",cGraph[i][0],deg);

		for(int j=2; j< deg+2;j++){
			fprintf(ogFile,":%d,%d",cGraph[i][j],1);
		}

		fprintf(ogFile,"\n");
	}

}

//argv[1]: input file
//argv[2]: number of vertices of input graph
//argv[3]: size of input graph
//argv[4]: partitioning method
int main(int argc, char * argv[]) {

	if(argc!=4){
		fprintf(stderr,"argv[1]-input file, argv[2]-number of vertices argv[3] size \n");
		exit(1);
	}

	iGraph = fopen(argv[1],"rb");
	if(iGraph == NULL){
		fprintf(stderr,"cannot open input file %s \n",argv[1]);
		exit(1);
	}

	Gno = atoi(argv[2]);
	Gsz = atol(argv[3]);
	//Gsz = 4409148160;

	rbuff = (char *)malloc(BLK_SZ);
	wbuff = (char *)malloc(BLK_SZ);
	mem = (char *)malloc((maxBlkPart+4)*BLK_SZ);
	cGraph = (vid **)malloc(maxNoPart*SZ_PTR);
	result = (count_type *)malloc(Gno*sizeof(count_type));
	adjBuff = (vid *)malloc((maxDeg+2)*SZ_VID);
	cmap = (vid *)malloc(Gno*SZ_VID);

	totalNumber = 0;
	noRead = 0;
	noWrite = 0;

	partTime = 0;

	Runtimecounter rt;
	rt.start();
	printf("maxBlkPart %d\n",maxBlkPart);
	emtcC();
	rt.stop();
	printf("%TIME total time: %f s \n",rt.GetRuntime());
	printf("partition Time %f s\n",partTime);
	printf("No. of fread %ld \n",noRead);
	printf("No. of fwrite %ld \n",noWrite);
	output("BTCTriA.txt");


	free(rbuff);
	free(wbuff);
	free(mem);
	free(cGraph);
	free(result);
	free(adjBuff);
	free(cmap);

	return 0;
}
