/*-------------------------------------------------------------------------
 *
 * CUDA function for backrpojection using FDK weigts for CBCT
 *
 *
 * CODE by  Ander Biguri
 * ---------------------------------------------------------------------------
 * ---------------------------------------------------------------------------
 * Copyright (c) 2015, University of Bath and CERN- European Organization for
 * Nuclear Research
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * ---------------------------------------------------------------------------
 *
 * Contact: tigre.toolbox@gmail.com
 * Codes  : https://github.com/CERN/TIGRE
 * ---------------------------------------------------------------------------
 */

#define  PI_2 1.57079632679489661923
#include <algorithm>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "voxel_backprojection_motion.hpp"
#include "mex.h"
#include <math.h>

// https://stackoverflow.com/questions/16282136/is-there-a-cuda-equivalent-of-perror
#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                mexErrMsgIdAndTxt("CBCT:CUDA:Atb",cudaGetErrorString(__err));\
        } \
} while (0)
    
    
    
#define MAXTREADS 1024
    /*GEOMETRY DEFINITION
     *
     *                Detector plane, behind
     *            |-----------------------------|
     *            |                             |
     *            |                             |
     *            |                             |
     *            |                             |
     *            |      +--------+             |
     *            |     /        /|             |
     *   A Z      |    /        / |*D           |
     *   |        |   +--------+  |             |
     *   |        |   |        |  |             |
     *   |        |   |     *O |  +             |
     *   *--->y   |   |        | /              |
     *  /         |   |        |/               |
     * V X        |   +--------+                |
     *            |-----------------------------|
     *
     *           *S
     *
     *
     *
     *
     *
     **/
texture<float, cudaTextureType3D , cudaReadModeElementType> tex;
texture<float, cudaTextureType3D , cudaReadModeElementType> tex_img;
texture<float, cudaTextureType3D , cudaReadModeElementType> texmvfX;
texture<float, cudaTextureType3D , cudaReadModeElementType> texmvfY;
texture<float, cudaTextureType3D , cudaReadModeElementType> texmvfZ;
surface<void,  cudaSurfaceType3D> surf_img;

const int PROJ_PER_KERNEL = 32;  // Number of 2D projections to be analyzed by a single thread. This can be tweaked to see what works best. 32 was the optimal value in the paper by Zinsser and Keck.
const int VOXELS_PER_THREAD = 8;  // Number of voxels to be computed by s single thread. Can be tweaked to see what works best. 4 was the optimal value in the paper by Zinsser and Keck.

// We have PROJ_PER_KERNEL projections and we need 6 parameters for each projection:
//   deltaX, deltaY, deltaZ, xyzOrigin, offOrig, offDetec
// So we need to keep PROJ_PER_KERNEL*6 values in our deltas array FOR EACH CALL to our main kernel
// (they will be updated in the main loop before each kernel call).

__constant__ Point3D projParamsArrayDev[6*PROJ_PER_KERNEL];  // Dev means it is on device

// We also need a corresponding array on the host side to be filled before each kernel call, then copied to the device (array in constant memory above)
Point3D projParamsArrayHost[6*PROJ_PER_KERNEL];   // Host means it is host memory

// Now we also need to store sinAlpha and cosAlpha for each projection (two floats per projection)
__constant__ float projSinCosArrayDev[3*PROJ_PER_KERNEL];

float projSinCosArrayHost[3*PROJ_PER_KERNEL];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// END RB, 10/31/2016: Add constant memory arrays to store parameters for all projections to be analyzed during a single kernel call
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//______________________________________________________________________________
//
//      Function:       rollPitchYaw
//
//      Description:    Main FDK backprojection kernel
//______________________________________________________________________________



//______________________________________________________________________________
//
//      Function:       kernelPixelBackprojectionFDK
//
//      Description:    Main FDK backprojection kernel
//______________________________________________________________________________

__global__ void kernelPixelBackprojectionFDK(const Geometry geo, float* image,const int currProjSetNumber, const int totalNoOfProjections,mfvInfo mvfData)
{
    
    // Old kernel call signature:
    // kernelPixelBackprojectionFDK<<<grid,block>>>(geo,dimage,i,deltaX,deltaY,deltaZ,xyzOrigin,offOrig,offDetec,sinalpha,cosalpha);
    // We just read in most of the params from the constant memory instead of getting them from the param list.
    // This is because we now have MANY params, since single kernel processes more than one projection!
    /* __global__ void kernelPixelBackprojectionFDK(const Geometry geo,
     * float* image,
     * const int indAlpha,
     * const Point3D deltaX ,
     * const Point3D deltaY,
     * const Point3D deltaZ,
     * const Point3D xyzOrigin,
     * const Point3D xyzOffset,
     * const Point3D uv0Offset,
     * const float sinalpha,
     * const float cosalpha){
     */
    unsigned long indY = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long indX = blockIdx.x * blockDim.x + threadIdx.x;
    // unsigned long startIndZ = blockIdx.z * blockDim.z + threadIdx.z;  // This is only STARTING z index of the column of voxels that the thread will handle
    unsigned long startIndZ = blockIdx.z * VOXELS_PER_THREAD + threadIdx.z;  // This is only STARTING z index of the column of voxels that the thread will handle
    //Make sure we dont go out of bounds
    if (indX>=geo.nVoxelX | indY>=geo.nVoxelY |startIndZ>=geo.nVoxelZ)
        return;
    
    
    
    float imageMVFratioX=(float)mvfData.nVoxelX/(float)geo.nVoxelX;
    float imageMVFratioY=(float)mvfData.nVoxelY/(float)geo.nVoxelY;
    float imageMVFratioZ=(float)mvfData.nVoxelZ/(float)geo.nVoxelZ;

    
    // We'll keep a local auxiliary array of values of a column of voxels that this thread will update
    float voxelColumn[VOXELS_PER_THREAD];
    
    // First we need to copy the curent 3D volume values from the column to our auxiliary array so that we can then
    // work on them (update them by computing values from multiple projections) locally - avoiding main memory reads/writes
    
    int colIdx;
    
    for(colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
    {
        unsigned long indZ = startIndZ + colIdx;
        // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
        // be trying to copy the out of bounds values back to the 3D volume anyway (bounds checks will be done in the final loop where the updated values go back to the main volume)
        if(indZ>=geo.nVoxelZ)
            break;   // break the loop.
        
        unsigned long long idx =indZ*geo.nVoxelX*geo.nVoxelY+indY*geo.nVoxelX + indX;
        voxelColumn[colIdx] = image[idx];   // Read the current volume value that we'll update by computing values from MULTIPLE projections (not just one)
        // We'll be updating the local (register) variable, avoiding reads/writes from the slow main memory.
    }  // END copy 3D volume voxels to local array
    
    // Now iterate through projections
    for(int projNumber=0; projNumber<PROJ_PER_KERNEL; projNumber++)
    {
        // Get the current parameters from parameter arrays in constant memory.
        int indAlpha = currProjSetNumber*PROJ_PER_KERNEL+projNumber;  // This is the ABSOLUTE projection number in the projection array
        
        // Our currImageVal will be updated by hovewer many projections we had left in the "remainder" - that's OK.
        if(indAlpha>=totalNoOfProjections)
            break;
        
        Point3D deltaX = projParamsArrayDev[6*projNumber];  // 6*projNumber because we have 6 Point3D values per projection
        Point3D deltaY = projParamsArrayDev[6*projNumber+1];
        Point3D deltaZ = projParamsArrayDev[6*projNumber+2];
        Point3D xyzOrigin = projParamsArrayDev[6*projNumber+3];
        Point3D xyzOffset = projParamsArrayDev[6*projNumber+4];
        Point3D S = projParamsArrayDev[6*projNumber+5];
        
        float sinalpha = projSinCosArrayDev[3*projNumber];     // 2*projNumber because we have 2 float (sin or cos angle) values per projection
        float cosalpha = projSinCosArrayDev[3*projNumber+1];
        float COR = projSinCosArrayDev[3*projNumber+2];
        
        // Geometric trasnformations:
        //Source, scaled XYZ coordinates

        // Now iterate through Z in our voxel column FOR A GIVEN PROJECTION
        for(colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
        {
            unsigned long indZ = startIndZ + colIdx;
            
            // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
            // be trying to copy the out of bounds values anyway (bounds checks will be done in the final loop where the values go to the main volume)
            if(indZ>=geo.nVoxelZ)
                break;   // break the loop.
            
                // "XYZ" in the scaled coordinate system of the current point. The image is rotated with the projection angles.
            Point3D P;
            P.x=(xyzOrigin.x+indX*deltaX.x+indY*deltaY.x+indZ*deltaZ.x);
            P.y=(xyzOrigin.y+indX*deltaX.y+indY*deltaY.y+indZ*deltaZ.y)-COR/geo.dDetecU;
            P.z=(xyzOrigin.z+indX*deltaX.z+indY*deltaY.z+indZ*deltaZ.z);
            
            // This is the vector defining the line from the source to the Voxel
            float vectX,vectY,vectZ;
            vectX=(P.x -S.x);
            vectY=(P.y -S.y);
            vectZ=(P.z -S.z);
            
            // Get the coordinates in the detector UV where the mid point of the voxel is projected.
            float t=(geo.DSO-geo.DSD /*-DOD*/ - S.x)/vectX;
            float y,z;
            y=vectY*t+S.y;
            z=vectZ*t+S.z;
            float u,v;
            u=y+geo.nDetecU/2-0.5;
            v=z+geo.nDetecV/2-0.5;
            
            float weigth;
            float realx,realy;
            realx=-geo.sVoxelX/2+geo.dVoxelX/2    +indX*geo.dVoxelX   +xyzOffset.x;
            realy=-geo.sVoxelY/2+geo.dVoxelY/2    +indY*geo.dVoxelY   +xyzOffset.y+COR;
            
            weigth=(geo.DSO+realy*sinalpha-realx*cosalpha)/geo.DSO;
            weigth=1/(weigth*weigth);
            
            // Get Value in the computed (U,V) and multiply by the corresponding weigth.
            // indAlpha is the ABSOLUTE number of projection in the projection array (NOT the current number of projection set!)
            voxelColumn[colIdx]+=tex3D(tex, u +0.5 ,
                    v +0.5 ,
                    indAlpha+0.5)*weigth;
            
        }  // END iterating through column of voxels
        
    }  // END iterating through multiple projections
    
    // And finally copy the updated local voxelColumn array back to our 3D volume (main memory)
    for(colIdx=0; colIdx<VOXELS_PER_THREAD; colIdx++)
    {
        unsigned long indZ = startIndZ + colIdx;
        // If we are out of bounds, break the loop. The voxelColumn array will be updated partially, but it is OK, because we won't
        // be trying to copy the out of bounds values back to the 3D volume anyway (bounds checks will be done in the final loop where the values go to the main volume)
        if(indZ>=geo.nVoxelZ)

            break;   // break the loop.
        
        unsigned long long idx =indZ*geo.nVoxelX*geo.nVoxelY+indY*geo.nVoxelX + indX;
        //image[idx] = voxelColumn[colIdx];   // Read the current volume value that we'll update by computing values from MULTIPLE projections (not just one)
        // We'll be updating the local (register) variable, avoiding reads/writes from the slow main memory.
        // According to references (Papenhausen), doing = is better than +=, since += requires main memory read followed by a write.
        // We did all the reads into the local array at the BEGINNING of this kernel. According to Papenhausen, this type of read-write split is
        // better for avoiding memory congestion.
        surf3Dwrite(voxelColumn[colIdx], surf_img, indX*sizeof(float), indY, indZ, cudaBoundaryModeZero);
    }  // END copy updated voxels from local array to our 3D volume
    
}  // END kernelPixelBackprojectionFDK


__global__ void kernelApplyMotion (const Geometry geo,
        float* image,
        mfvInfo mvfData){
    unsigned long indY = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned long indX = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long indZ = blockIdx.z * blockDim.z + threadIdx.z;
    //Make sure we dont go out of bounds
    unsigned long idx =indZ*geo.nVoxelX*geo.nVoxelY+indY*geo.nVoxelX + indX;
    if (indX>=geo.nVoxelX | indY>=geo.nVoxelY |indZ>=geo.nVoxelZ)
        return;
    
    
    float imageMVFratioX=(float)mvfData.nVoxelX/(float)geo.nVoxelX;
    float imageMVFratioY=(float)mvfData.nVoxelY/(float)geo.nVoxelY;
    float imageMVFratioZ=(float)mvfData.nVoxelZ/(float)geo.nVoxelZ;
    float xm=tex3D(texmvfX,(indX*imageMVFratioX+0.5)/mvfData.nVoxelX, (indY*imageMVFratioY+0.5)/mvfData.nVoxelY, (indZ*imageMVFratioZ+0.5)/mvfData.nVoxelZ);
    float ym=tex3D(texmvfY,(indX*imageMVFratioX+0.5)/mvfData.nVoxelX, (indY*imageMVFratioY+0.5)/mvfData.nVoxelY, (indZ*imageMVFratioZ+0.5)/mvfData.nVoxelZ);
    float zm=tex3D(texmvfZ,(indX*imageMVFratioX+0.5)/mvfData.nVoxelX, (indY*imageMVFratioY+0.5)/mvfData.nVoxelY, (indZ*imageMVFratioZ+0.5)/mvfData.nVoxelZ);
    image[idx]+=tex3D(tex_img,indX+xm/geo.dVoxelX+0.5,indY+ym/geo.dVoxelY+0.5,indZ+zm/geo.dVoxelZ+0.5);
    //image[idx]=zm;
}

int voxel_backprojection_motion(float const * const projections, Geometry geo, float* result,float const * const alphas,int nalpha,float const * const mvfX, float const * const mvfY,
        float const * const mvfZ, mfvInfo mvfData){
    
    
    /*
     * Allocate texture memory on the device
     */
    // copy data to CUDA memory
    cudaArray *d_projectiondata = 0;
    const cudaExtent extent = make_cudaExtent(geo.nDetecU,geo.nDetecV,nalpha);
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaMalloc3DArray(&d_projectiondata, &channelDesc, extent);
    cudaCheckErrors("cudaMalloc3D error 3D tex");

    cudaMemcpy3DParms copyParams = { 0 };
    copyParams.srcPtr = make_cudaPitchedPtr((void*)projections, extent.width*sizeof(float), extent.width, extent.height);
    copyParams.dstArray = d_projectiondata;
    copyParams.extent = extent;
    copyParams.kind = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);

    cudaCheckErrors("cudaMemcpy3D fail");
    
    // Configure texture options
    tex.normalized = false;
    tex.filterMode = cudaFilterModeLinear;
    tex.addressMode[0] = cudaAddressModeBorder;
    tex.addressMode[1] = cudaAddressModeBorder;
    tex.addressMode[2] = cudaAddressModeBorder;
    
    cudaBindTextureToArray(tex, d_projectiondata, channelDesc);
    
    cudaCheckErrors("3D texture memory bind fail");
            

    // copy MVF
    
    
    cudaArray *d_mvfX = 0;
    cudaArray *d_mvfY = 0;
    cudaArray *d_mvfZ = 0;
    
    const cudaExtent extent_mvf = make_cudaExtent(mvfData.nVoxelX, mvfData.nVoxelY, mvfData.nVoxelZ);
    cudaChannelFormatDesc channelDescX = cudaCreateChannelDesc<float>();
    cudaChannelFormatDesc channelDescY = cudaCreateChannelDesc<float>();
    cudaChannelFormatDesc channelDescZ = cudaCreateChannelDesc<float>();
    
    cudaMalloc3DArray(&d_mvfX, &channelDescX, extent_mvf);
    cudaMalloc3DArray(&d_mvfY, &channelDescY, extent_mvf);
    cudaMalloc3DArray(&d_mvfZ, &channelDescZ, extent_mvf);
    
    cudaCheckErrors("cudaMalloc3D error 3D tex");
    copyParams.srcPtr = make_cudaPitchedPtr((void*)mvfX, extent_mvf.width*sizeof(float), extent_mvf.width, extent_mvf.height);
    copyParams.dstArray = d_mvfX;
    copyParams.extent = extent_mvf;
    copyParams.kind = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);
    //
    copyParams.srcPtr = make_cudaPitchedPtr((void*)mvfY, extent_mvf.width*sizeof(float), extent_mvf.width, extent_mvf.height);
    copyParams.dstArray = d_mvfY;
    copyParams.extent = extent_mvf;
    copyParams.kind = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);
    //
    copyParams.srcPtr = make_cudaPitchedPtr((void*)mvfZ, extent_mvf.width*sizeof(float), extent_mvf.width, extent_mvf.height);
    copyParams.dstArray = d_mvfZ;
    copyParams.extent = extent_mvf;
    copyParams.kind = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);
    cudaCheckErrors("cudaMemcpy3D fail");
    // Configure texture options
    texmvfX.normalized = true;                           texmvfY.normalized = true;                       texmvfZ.normalized = true;
    texmvfX.filterMode = cudaFilterModeLinear;           texmvfY.filterMode = cudaFilterModeLinear;       texmvfZ.filterMode = cudaFilterModeLinear;
    texmvfX.addressMode[0] = cudaAddressModeBorder;      texmvfY.addressMode[0] = cudaAddressModeBorder;  texmvfZ.addressMode[0] = cudaAddressModeBorder;
    texmvfX.addressMode[1] = cudaAddressModeBorder;      texmvfY.addressMode[1] = cudaAddressModeBorder;  texmvfZ.addressMode[1] = cudaAddressModeBorder;
    texmvfX.addressMode[2] = cudaAddressModeBorder;      texmvfY.addressMode[2] = cudaAddressModeBorder;  texmvfZ.addressMode[2] = cudaAddressModeBorder;
    
    
    cudaBindTextureToArray(texmvfX, d_mvfX, channelDescX);
    cudaBindTextureToArray(texmvfY, d_mvfY, channelDescY);
    cudaBindTextureToArray(texmvfZ, d_mvfZ, channelDescZ);
    
    cudaCheckErrors("3D texture memory bind fail");
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Image and auxiliary image memory. Image is output, auxiliary is surf&tex.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Out image
    float *d_image = 0;
    size_t num_bytes = geo.nVoxelX * geo.nVoxelY * geo.nVoxelZ * sizeof(float);
    cudaMalloc((void**)&d_image, num_bytes);
    cudaMemset(d_image, 0, num_bytes);
    cudaCheckErrors("malloc d_image");
    // auxiliar
    
    cudaArray* d_image_aux;
    const cudaExtent extent_img = make_cudaExtent(geo.nVoxelX,  geo.nVoxelY, geo.nVoxelZ);
    cudaChannelFormatDesc channelDesc_img = cudaCreateChannelDesc<float>();
    cudaMalloc3DArray(&d_image_aux, &channelDesc, extent_img, cudaArraySurfaceLoadStore);

    tex_img.normalized = false;
    tex_img.filterMode = cudaFilterModeLinear;
    tex_img.addressMode[0] = cudaAddressModeBorder;
    tex_img.addressMode[1] = cudaAddressModeBorder;
    tex_img.addressMode[2] = cudaAddressModeBorder;
    // Bind image_aux to both texture and surface
    cudaBindSurfaceToArray(surf_img, d_image_aux, channelDesc_img);
    cudaBindTextureToArray(tex_img, d_image_aux, channelDesc_img);
    
    
    
    // If we are going to time
    bool timekernel=false;
    cudaEvent_t start, stop;
    float elapsedTime;
    
    
    int divx,divy,divz;
    int divxm,divym,divzm;
    // RB: Use the optimal (in their tests) block size from paper by Zinsser and Keck (16 in x and 32 in y).
    // I tried different sizes and shapes of blocks (tiles), but it does not appear to significantly affect trhoughput, so
    // let's stick with the values from Zinsser and Keck.
    divx=16;
    divy=32;
    divz=VOXELS_PER_THREAD;      // We now only have 32 x 16 threads per block (flat tile, see below), BUT each thread works on a Z column of VOXELS_PER_THREAD voxels, so we effectively need fewer blocks!
    dim3 grid((geo.nVoxelX+divx-1)/divx,
            (geo.nVoxelY+divy-1)/divy,
            (geo.nVoxelZ+divz-1)/divz);
    
    dim3 block(divx,divy,1);    // Note that we have 1 in the Z size, not divz, since each thread works on a vertical set of VOXELS_PER_THREAD voxels (so we only need a "flat" tile of threads, with depth of 1)
    
    divxm=32;
    divym=32;
    divzm=1;
    dim3 gridm((geo.nVoxelX+divxm-1)/divxm,
            (geo.nVoxelY+divym-1)/divym,
            (geo.nVoxelZ+divzm-1)/divzm);
    dim3 blockm(divxm,divym,divzm);
    //////////////////////////////////////////////////////////////////////////////////////
    // Main reconstruction loop: go through projections (rotation angles) and backproject
    //////////////////////////////////////////////////////////////////////////////////////
    
    // Since we'll have multiple projections processed by a SINGLE kernel call, compute how many
    // kernel calls we'll need altogether.
    int noOfKernelCalls = (nalpha+PROJ_PER_KERNEL-1)/PROJ_PER_KERNEL;  // We'll take care of bounds checking inside the loop if nalpha is not divisible by PROJ_PER_KERNEL
    
    for (unsigned int i=0; i<noOfKernelCalls; i++)
    {
        // Now we need to generate and copy all data for PROJ_PER_KERNEL projections to constant memory so that our kernel can use it
        int j;
        for(j=0; j<PROJ_PER_KERNEL; j++)
        {
            int currProjNumber=i*PROJ_PER_KERNEL+j;
            
            if(currProjNumber>=nalpha)
                break;  // Exit the loop. Even when we leave the param arrays only partially filled, this is OK, since the kernel will check bounds anyway.
            
            Point3D deltaX,deltaY,deltaZ,xyzOrigin, offOrig, /*offDetec,*/source;
            float sinalpha,cosalpha;
            
            geo.alpha=-alphas[currProjNumber];
            sinalpha=sin(geo.alpha);
            cosalpha=cos(geo.alpha);
            
            projSinCosArrayHost[3*j]=sinalpha;  // 2*j because we have 2 float (sin or cos angle) values per projection
            projSinCosArrayHost[3*j+1]=cosalpha;
            projSinCosArrayHost[3*j+2]=geo.COR;
            
            computeDeltasCube(geo,geo.alpha,currProjNumber,&xyzOrigin,&deltaX,&deltaY,&deltaZ,&source);
            
            offOrig.x=geo.offOrigX[currProjNumber];
            offOrig.y=geo.offOrigY[currProjNumber];
            
            
            projParamsArrayHost[6*j]=deltaX;		// 6*j because we have 6 Point3D values per projection
            projParamsArrayHost[6*j+1]=deltaY;
            projParamsArrayHost[6*j+2]=deltaZ;
            projParamsArrayHost[6*j+3]=xyzOrigin;
            projParamsArrayHost[6*j+4]=offOrig;
            projParamsArrayHost[6*j+5]=source;
        }   // END for (preparing params for kernel call)
        
        // Copy the prepared parameter arrays to constant memory to make it available for the kernel
        cudaMemcpyToSymbol(projSinCosArrayDev, projSinCosArrayHost, sizeof(float)*3*PROJ_PER_KERNEL);
        cudaMemcpyToSymbol(projParamsArrayDev, projParamsArrayHost, sizeof(Point3D)*6*PROJ_PER_KERNEL);
        if (timekernel){
            cudaEventCreate(&start);
            cudaEventRecord(start,0);
        }
        kernelPixelBackprojectionFDK<<<grid,block>>>(geo,d_image,i,nalpha,mvfData);
        cudaCheckErrors("Kernel fail 1");
       if (timekernel){
            cudaEventCreate(&stop);
            cudaEventRecord(stop,0);
            cudaEventSynchronize(stop);
            cudaEventElapsedTime(&elapsedTime, start,stop);
            mexPrintf("%f\n" ,elapsedTime);
            cudaCheckErrors("cuda Timing fail");
            
        }
        kernelApplyMotion<<<gridm,blockm>>>(geo, d_image, mvfData);
        cudaCheckErrors("Kernel fail 2");
        
         
    }  
        
    
    cudaMemcpy(result, d_image, num_bytes, cudaMemcpyDeviceToHost);
    cudaCheckErrors("cudaMemcpy result fail");
    
    cudaUnbindTexture(tex);
    cudaUnbindTexture(texmvfX);
    cudaUnbindTexture(texmvfY);
    cudaUnbindTexture(texmvfZ);
    cudaUnbindTexture(tex_img);
    
    
    cudaCheckErrors("Unbind  fail");
    
    cudaFree(d_image);
    cudaFreeArray(d_mvfX);
    cudaFreeArray(d_mvfY);
    cudaFreeArray(d_mvfZ);
    
    cudaFreeArray(d_projectiondata);
    cudaFreeArray(d_image_aux);
    cudaCheckErrors("cudaFree d_imagedata fail");
    //cudaDeviceReset();

    return 0;
    
}
void computeDeltasCube(Geometry geo, float alpha,int i, Point3D* xyzorigin, Point3D* deltaX, Point3D* deltaY, Point3D* deltaZ,Point3D* S)
{
    Point3D P0, Px0,Py0,Pz0, source;
    // Get coords of Img(0,0,0)
    P0.x=-(geo.sVoxelX/2-geo.dVoxelX/2)+geo.offOrigX[i];
    P0.y=-(geo.sVoxelY/2-geo.dVoxelY/2)+geo.offOrigY[i];
    P0.z=-(geo.sVoxelZ/2-geo.dVoxelZ/2)+geo.offOrigZ[i];
    
    // Get coors from next voxel in each direction
    Px0.x=P0.x+geo.dVoxelX;       Py0.x=P0.x;                Pz0.x=P0.x;
    Px0.y=P0.y;                   Py0.y=P0.y+geo.dVoxelY;    Pz0.y=P0.y;
    Px0.z=P0.z;                   Py0.z=P0.z;                Pz0.z=P0.z+geo.dVoxelZ;
    
    // Rotate image (this is equivalent of rotating the source and detector)
    
    Point3D P, Px,Py,Pz; // We need other auxiliar variables to be able to perform the rotation, or we would overwrite values!
    P.x =P0.x *cos(alpha)-P0.y *sin(alpha);       P.y =P0.x *sin(alpha)+P0.y *cos(alpha);      P.z =P0.z;
    Px.x=Px0.x*cos(alpha)-Px0.y*sin(alpha);       Px.y=Px0.x*sin(alpha)+Px0.y*cos(alpha);      Px.z=Px0.z;
    Py.x=Py0.x*cos(alpha)-Py0.y*sin(alpha);       Py.y=Py0.x*sin(alpha)+Py0.y*cos(alpha);      Py.z=Py0.z;
    Pz.x=Pz0.x*cos(alpha)-Pz0.y*sin(alpha);       Pz.y=Pz0.x*sin(alpha)+Pz0.y*cos(alpha);      Pz.z=Pz0.z;
    
    //detector offset
    P.z =P.z-geo.offDetecV[i];          P.y =P.y-geo.offDetecU[i];
    Px.z =Px.z-geo.offDetecV[i];          Px.y =Px.y-geo.offDetecU[i];
    Py.z =Py.z-geo.offDetecV[i];          Py.y =Py.y-geo.offDetecU[i];
    Pz.z =Pz.z-geo.offDetecV[i];          Pz.y =Pz.y-geo.offDetecU[i];
    //Detector Roll pitch Yaw
    //
    //
    // first, we need to offset everything so (0,0,0) is the center of the detector
    // Only X is required for that
    P.x=P.x+(geo.DSD-geo.DSO);
    Px.x=Px.x+(geo.DSD-geo.DSO);
    Py.x=Py.x+(geo.DSD-geo.DSO);
    Pz.x=Pz.x+(geo.DSD-geo.DSO);
    rollPitchYawT(geo,i,&P);
    rollPitchYawT(geo,i,&Px);
    rollPitchYawT(geo,i,&Py);
    rollPitchYawT(geo,i,&Pz);
    
    P.x=P.x-(geo.DSD-geo.DSO);
    Px.x=Px.x-(geo.DSD-geo.DSO);
    Py.x=Py.x-(geo.DSD-geo.DSO);
    Pz.x=Pz.x-(geo.DSD-geo.DSO);
    //Done for P, now source
    
    source.x=geo.DSD; //allready offseted for rotation of teh detector
    source.y=-geo.offDetecU[i];
    source.z=-geo.offDetecV[i];
    rollPitchYawT(geo,i,&source);
    
    source.x=source.x-(geo.DSD-geo.DSO);
    
    // Scale coords so detector pixels are 1x1
    
    P.z =P.z /geo.dDetecV;                          P.y =P.y/geo.dDetecU;
    Px.z=Px.z/geo.dDetecV;                          Px.y=Px.y/geo.dDetecU;
    Py.z=Py.z/geo.dDetecV;                          Py.y=Py.y/geo.dDetecU;
    Pz.z=Pz.z/geo.dDetecV;                          Pz.y=Pz.y/geo.dDetecU;
    source.z=source.z/geo.dDetecV;                  source.y=source.y/geo.dDetecU;
    // get deltas of the changes in voxels
    deltaX->x=Px.x-P.x;   deltaX->y=Px.y-P.y;    deltaX->z=Px.z-P.z;
    deltaY->x=Py.x-P.x;   deltaY->y=Py.y-P.y;    deltaY->z=Py.z-P.z;
    deltaZ->x=Pz.x-P.x;   deltaZ->y=Pz.y-P.y;    deltaZ->z=Pz.z-P.z;
    
    
    *xyzorigin=P;
    *S=source;
}  // END computeDeltasCube
void rollPitchYawT(Geometry geo,int i, Point3D* point){
    Point3D auxPoint;
    auxPoint.x=point->x;
    auxPoint.y=point->y;
    auxPoint.z=point->z;
    
    point->x=cos(geo.dRoll[i])*cos(geo.dPitch[i])*auxPoint.x
            +sin(geo.dRoll[i])*cos(geo.dPitch[i])*auxPoint.y
            -sin(geo.dPitch[i])*auxPoint.z;
    
    
    point->y=(cos(geo.dRoll[i])*sin(geo.dPitch[i])*sin(geo.dYaw[i]) - sin(geo.dRoll[i])*cos(geo.dYaw[i]))*auxPoint.x
            +(sin(geo.dRoll[i])*sin(geo.dPitch[i])*sin(geo.dYaw[i]) + cos(geo.dRoll[i])*cos(geo.dYaw[i]))*auxPoint.y
            +cos(geo.dPitch[i])*sin(geo.dYaw[i])*auxPoint.z;
    
    
    point->z=(cos(geo.dRoll[i])*sin(geo.dPitch[i])*cos(geo.dYaw[i]) + sin(geo.dRoll[i])*sin(geo.dYaw[i]))*auxPoint.x
            +(sin(geo.dRoll[i])*sin(geo.dPitch[i])*cos(geo.dYaw[i]) - cos(geo.dRoll[i])*sin(geo.dYaw[i]))*auxPoint.y
            +cos(geo.dPitch[i])*cos(geo.dYaw[i])*auxPoint.z;
    
}
