#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "com_falconcomputing_genomics_haplotypecaller_FalconPairhmm.h"

#include "pairhmm/common/pairhmm_common.h"
#include "pairhmm/avx/avx_impl.h"
#include "pairhmm/common/Context.h"
#include "pairhmm/client/PairhmmClient.h"
#include "pairhmm/common/PairhmmJavaData.h"
#include <xmmintrin.h>



JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_initNative
(JNIEnv* env, jclass cls, jclass readDataHolder, jclass haplotypeDataHolder, jboolean use_double, jint max_threads, jboolean use_fpga){
    struct PairhmmContext* curCtx = init_pairhmm_context();
    JavaData javaData;
    javaData.init(env, readDataHolder, haplotypeDataHolder);
    curCtx->g_use_double = use_double;
    curCtx->g_use_fpga = use_fpga;

    // enable FTZ
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

    // init convert char table
    ConvertChar::init();
    for(int i = 0; i < 3; i++)
        inputBank[i].first = "bankID";
#ifdef DEPLOY_aws
    inputBank[0].second = 3; inputBank[1].second = 1; inputBank[2].second = 2;
#else
    inputBank[0].second = 0; inputBank[1].second = 1; inputBank[2].second = 3;
#endif
}

void compute_violate(pairhmmInput* input, double* results, bool use_double){
    struct PairhmmContext* curCtx = init_pairhmm_context();
    struct timespec time1, time2, diff;
    clock_gettime(CLOCK_REALTIME, &time1);
    testcase testCase;
    int index = 0;
    for (int r = 0; (size_t)r < input->reads.size(); r++){
        testCase.rslen = input->reads[r].bases.size(); 
        testCase.rs = input->reads[r].bases.c_str(); 
        testCase.i = input->reads[r]._i.c_str();
        testCase.d = input->reads[r]._d.c_str();
        testCase.c = input->reads[r]._c.c_str();
        testCase.q = input->reads[r]._q.c_str();
        for(int j = 0; (size_t)j < input->haps.size(); j++){
            testCase.haplen = input->haps[j].bases.size();
            testCase.hap = input->haps[j].bases.c_str();
            double result_final = 0;
            float result_float = use_double ? 0.0f : compute_fp_avxs(&testCase);
            if(result_float < MIN_ACCEPTED){
                double result_double = compute_fp_avxd(&testCase);
                result_final = log10(result_double) - g_ctxd.LOG10_INITIAL_CONSTANT;
            }
            else{
                result_final = (double)(log10f(result_float) - g_ctxf.LOG10_INITIAL_CONSTANT);
            }
            results[index++] = result_final;
        }
    }
    clock_gettime(CLOCK_REALTIME, &time2);
    diff = pmm_diff_time(time1, time2);
    curCtx->direct_avx_time += diff.tv_sec * 1e9 + diff.tv_nsec;

}

double countCell(pairhmmInput* input, short maxCols, bool& violate){    
    double numCell = 0;
    for(int i = 0; (size_t)i < input->reads.size(); i++){
        for(int j = 0; (size_t)j < input->haps.size(); j++){
            int cur_numCell = input->reads[i].bases.size() * input->haps[j].bases.size();
            if(input->reads[i].bases.size() * (input->reads[i].bases.size() + maxCols) <= 1 + READ_BLOCK_SIZE * MAX_READ_LEN){
                violate = true;
                return numCell;
            }
            numCell += cur_numCell;
        }
    }
    return numCell;
}

bool worthFPGA(pairhmmInput* input, short maxCols, double cellNum){
    if(input->reads.size() < DIE_NUM)
        return false;
    bool violate = false;
    if(cellNum <= 0)
        cellNum = countCell(input, maxCols, violate);
    if(violate)
        return false;
    double AVX_time = cellNum / AVX_perf; //unit nano seconds, 0.6 GCUPS
    double CPU_to_DRAM = DIE_NUM * 300000 + DIE_NUM * sizeof(FPGAInput) / 3.0;
    double DRAM_to_FPGA = 0;
    double DRAM_to_CPU = DIE_NUM * 300000 + (input->reads.size() * input->haps.size() * 4) / 3.0;
    double FPGA_compute = cellNum / FPGA_perf;
    if(AVX_time > (CPU_to_DRAM + DRAM_to_FPGA + DRAM_to_CPU + FPGA_compute)){
        return true;
    }
    else
        return false;
}

void convert_read_input(readDataPack* cur_host_read, Read* input, int idx){
    uint8_t _rs = ConvertChar::get(input->bases[idx]);
    uint8_t _q = input->_q[idx];
    uint8_t _i = input->_i[idx];
    uint8_t _d = input->_d[idx];
    uint8_t _c = input->_c[idx];
    uint32_t score = ((uint32_t)(_rs & 0x7) << 28) | ((uint32_t)(_q & 127) << 21) | \
                     ((uint32_t)(_i & 127) << 14) | ((uint32_t)(_d & 127) << 7) | \
                     (uint32_t)(_c & 127);
    float m2m_val = m2m_table[128 * _i + _d];
    cur_host_read->scores = score;
    cur_host_read->m2m = m2m_val;
}

void distributeReads(pairhmmInput* input, int read_base_index, int hap_base_index, int& numRead0, int& numRead1, int& numRead2, int& totalNumRead, int& totalNumHap, short maxCols, bool& violate){
    //first get the total actual cells
    double numCell = 0;
    for(int i = read_base_index; i < read_base_index + totalNumRead; i++){
        for(int j = hap_base_index; j < hap_base_index + totalNumHap; j++){
            int amended_read_length = input->reads[i].bases.size() + 1;
            int new_rows = amended_read_length;
            if(new_rows < DEP_DIST)
                new_rows = DEP_DIST;
            int cur_numCell = (new_rows + 1) * (amended_read_length + maxCols);
            numCell += cur_numCell;
        }
    }
    if(!worthFPGA(input, maxCols, numCell)){
        violate = true;
        return;
    }
    float SLR_numCells[3];
    if(DIE_NUM == 2){
        SLR_numCells[0] = floor(numCell * SLR0_PE_NUM / (SLR0_PE_NUM + SLR1_PE_NUM));
        SLR_numCells[1] = numCell - SLR_numCells[0];
        SLR_numCells[2] = 0;
    }
    else{
        SLR_numCells[0] = floor(numCell * SLR0_PE_NUM / (SLR0_PE_NUM + SLR1_PE_NUM + SLR2_PE_NUM));
        SLR_numCells[1] = floor(numCell * SLR1_PE_NUM / (SLR0_PE_NUM + SLR1_PE_NUM + SLR2_PE_NUM));
        SLR_numCells[2] = numCell - SLR_numCells[0] - SLR_numCells[1];
    }
    float curCells = 0.0;
    int readCount = 0;
    int i = 0;
    for(i = read_base_index; i < read_base_index + totalNumRead; i++){
        for(int j = hap_base_index; j < hap_base_index + totalNumHap; j++){
            int amended_read_length = input->reads[i].bases.size() + 1;
            int new_rows = amended_read_length;
            if(new_rows < DEP_DIST)
                new_rows = DEP_DIST;
            curCells += (new_rows + 1) * (amended_read_length + maxCols);
        }
        readCount++;
        if(curCells >= SLR_numCells[0]){
            numRead0 = readCount;
            break;
        }
    }
    if(numRead0 <= 0){
        violate = true;
        return;
    }
    if(DIE_NUM == 2){
        numRead1 = totalNumRead - numRead0;
        if(numRead1 <= 0)
            violate = true;
        numRead2 = 0;
        return;
    }
    else{
        curCells = 0;
        readCount = 0;
        i++;
        for( ; i < read_base_index + totalNumRead; i++){
            for(int j = hap_base_index; j < hap_base_index + totalNumHap; j++){
                int amended_read_length = input->reads[i].bases.size() + 1;
                int new_rows = amended_read_length;
                if(new_rows < DEP_DIST)
                    new_rows = DEP_DIST;
                curCells += (new_rows + 1) * (amended_read_length + maxCols);
            }
            readCount++;
            if(curCells >= SLR_numCells[1]){
                numRead1 = readCount;
                break;
            }
        }
        numRead2 = totalNumRead - numRead0 - numRead1;
        if(numRead2 <= 0)
            violate = true;
    }

}

void computePairhmmAVXSegment(pairhmmInput* input, int read_base_index, int hap_base_index, int cur_numRead, int cur_numHap, vector<float>& output){
    struct PairhmmContext* curCtx = init_pairhmm_context();
    struct timespec time1, time2, diff;
    clock_gettime(CLOCK_REALTIME, &time1);
    testcase testCase;
    for(int i = read_base_index; i < read_base_index + cur_numRead; i++){
        testCase.rslen = input->reads[i].bases.size();
        testCase.rs = input->reads[i].bases.c_str();
        testCase.i = input->reads[i]._i.c_str();
        testCase.d = input->reads[i]._d.c_str();
        testCase.c = input->reads[i]._c.c_str();
        testCase.q = input->reads[i]._q.c_str();
        for(int j = hap_base_index; j < hap_base_index + cur_numHap; j++){
            testCase.haplen = input->haps[j].bases.size();
            testCase.hap = input->haps[j].bases.c_str();
            float result_final = compute_fp_avxs(&testCase);
            output[i * input->haps.size() + j] = result_final;
        }
    }
    clock_gettime(CLOCK_REALTIME, &time2);
    diff = pmm_diff_time(time1, time2);
    curCtx->direct_avx_time += diff.tv_sec * 1e9 + diff.tv_nsec;


}

bool cmpReadInfo(struct readInfo a, struct readInfo b){
    return (a.new_rows > b.new_rows);
}

void sortReads(pairhmmInput* input, int read_base_index, int cur_numRead, int cur_numHap, short maxCols, int slr_pu_num[3], FPGAInput* host_input[3]){
    vector<readInfo> sortedReadInfo[3];
    int read_start_index = 0; 
    for(int k = 0; k < DIE_NUM; k++){
        for(int i = 0; i < host_input[k]->numRead; i += 2){
            readInfo curInfo;
            int cur_read_len = input->reads[i + read_start_index + read_base_index].bases.size();
            curInfo.readLen[0] = cur_read_len;
            curInfo.big_rows = cur_read_len + 1;
            curInfo.oneOrTwo = false;
            curInfo.resultOffset = cur_numHap * i;
            curInfo.readID[0] = i;
            if(i + 1 < host_input[k]->numRead){
                cur_read_len = input->reads[i + 1 + read_start_index + read_base_index].bases.size();
                curInfo.readLen[1] = cur_read_len;
                curInfo.oneOrTwo = true;
                if(cur_read_len + 1 > curInfo.big_rows)
                    curInfo.big_rows = cur_read_len + 1;
                curInfo.readID[1] = i + 1;
            }
            if(curInfo.big_rows < DEP_DIST)
                curInfo.new_rows = DEP_DIST;
            else
                curInfo.new_rows = curInfo.big_rows;
            curInfo.curIterNum = (curInfo.new_rows + 1) * (curInfo.big_rows + maxCols);
            curInfo.infoPacked = curInfo.readLen[0] + (curInfo.readLen[1] << 8) + ((uint64_t)curInfo.resultOffset << 16) + (((uint64_t)curInfo.curIterNum - 1) << 37) + ((uint64_t)curInfo.oneOrTwo << 58);
            sortedReadInfo[k].push_back(curInfo);
        }
        //printf("sortedReadInfo[%d] has size %d\n", k, sortedReadInfo[k].size());
        sort(sortedReadInfo[k].begin(), sortedReadInfo[k].end(), cmpReadInfo);
        for(int j = 0; j < sortedReadInfo[k].size(); j+=slr_pu_num[k]){
            int upper_bound = slr_pu_num[k];
            if(j + upper_bound >= sortedReadInfo[k].size())
                upper_bound = sortedReadInfo[k].size() - j;
            for(int m = 0; m < upper_bound / 2; m++){
                swap(sortedReadInfo[k][j + m], sortedReadInfo[k][j + upper_bound - 1 - m]); 
            }
        }
        read_start_index += host_input[k]->numRead;
    }
    read_start_index = 0;
    for(int slr_id = 0; slr_id < 3; slr_id++){
        int PU_id = 0;
        for(int i = 0; i < sortedReadInfo[slr_id].size(); i++){
            for(int j = 0; j < sortedReadInfo[slr_id][i].readLen[0]; j++){
                convert_read_input(&(host_input[slr_id]->dataPack.readData[READ_BLOCK_SIZE * i][j]), &(input->reads[sortedReadInfo[slr_id][i].readID[0] + read_base_index + read_start_index]), j);
            }
            if(sortedReadInfo[slr_id][i].oneOrTwo){
                for(int j = 0; j < sortedReadInfo[slr_id][i].readLen[1]; j++){
                    convert_read_input(&(host_input[slr_id]->dataPack.readData[READ_BLOCK_SIZE * i + 1][j]), &(input->reads[sortedReadInfo[slr_id][i].readID[1] + read_base_index + read_start_index]), j);
                }
            }

        }
        read_start_index += host_input[slr_id]->numRead;
    }

    for(int slr_id = 0; slr_id < 3; slr_id++){
        int PU_id = 0;
        for(int i = 0; i < sortedReadInfo[slr_id].size(); i++){
            host_input[slr_id]->dataPack.iterNum[PU_id] += sortedReadInfo[slr_id][i].curIterNum;
            host_input[slr_id]->dataPack.readInfo[i] = sortedReadInfo[slr_id][i].infoPacked;
            host_input[slr_id]->dataPack.numReadPU[PU_id] += (sortedReadInfo[slr_id][i].oneOrTwo) ? 2 : 1;
            PU_id = (PU_id == slr_pu_num[slr_id] - 1) ? 0 : PU_id + 1;

        }
    }

    int hap_batch_size = cur_numHap / HAP_BLOCK_SIZE + (cur_numHap % HAP_BLOCK_SIZE != 0);
    for(int slr_id = 0; slr_id < 3; slr_id++){
        long long total_slr_itercount = 0;
        for(int k = 0; k < slr_pu_num[slr_id]; k++){
            host_input[slr_id]->dataPack.iterNum[k] *= hap_batch_size;
            total_slr_itercount += host_input[slr_id]->dataPack.iterNum[k];
        }
    }
    
    
}


float update_host_inputs(pairhmmInput* input, int read_base_index, int hap_base_index, int cur_numRead, int cur_numHap, bool& violate, FPGAInput* host_input[3]){
    float cells = 0.0;
    for(int i = read_base_index; i < cur_numRead + read_base_index; i++){
        for(int j = hap_base_index; j < cur_numHap + hap_base_index; j++){
            cells = cells + input->reads[i].bases.size() * input->haps[j].bases.size();
        }
    }
    for(int i = 0; i < DIE_NUM; i++){
        host_input[i]->numHap = cur_numHap;
    }

    short maxCols = 0;
    for(int i = 0; i < cur_numHap; i++){
        if(input->haps[i + hap_base_index].bases.size() > MAX_HAP_LEN){
            violate = true;
            return 0;
        }
        short hapLenTmp = input->haps[i + hap_base_index].bases.size(); 
        if(hapLenTmp + 1 > maxCols){
            maxCols = hapLenTmp + 1;
        }
        for(int j = 0; j < DIE_NUM; j++){
            host_input[j]->dataPack.hapDataLen[i].hapLen = hapLenTmp;
            host_input[j]->dataPack.hapDataLen[i].oneDivHapLen = g_ctxf.INITIAL_CONSTANT / (float)(hapLenTmp);
        }
    }

    for(int i = 0; i < MAX_HAPDATA_NUM / HAP_BLOCK_SIZE; i++){
        for(int j = 0; j < MAX_HAP_LEN; j++){
            uint32_t tmp = 0;
            for(int k = 0; k < HAP_BLOCK_SIZE; k++){
                char cur_base;
                if(i * HAP_BLOCK_SIZE + k + hap_base_index >= (int)input->haps.size())
                    cur_base = 'N';
                else if((size_t)j >= input->haps[i * HAP_BLOCK_SIZE + k + hap_base_index].bases.size())
                    cur_base = 'N';
                else{
                    cur_base = input->haps[i * HAP_BLOCK_SIZE + k + hap_base_index].bases[j];
                }
                tmp = tmp | ((uint16_t)(ConvertChar::get(cur_base)) << (4 * k));
            }
            for(int k = 0; k < DIE_NUM; k++){
                host_input[k]->dataPack.hapData[i][j] = tmp;
            }
        }
    }

    //distribute reads to each PE
    int slr_pu_num[3];
    slr_pu_num[0] = SLR0_PE_NUM / (READ_BLOCK_SIZE * HAP_BLOCK_SIZE);
    slr_pu_num[1] = SLR1_PE_NUM / (READ_BLOCK_SIZE * HAP_BLOCK_SIZE);
    slr_pu_num[2] = SLR2_PE_NUM / (READ_BLOCK_SIZE * HAP_BLOCK_SIZE);
    for(int i = 0; i < DIE_NUM; i++){
        for(int j = 0; j < 64; j++){
            host_input[i]->dataPack.numReadPU[j] = 0;
            host_input[i]->dataPack.iterNum[j] = 0;
        }
    }

    distributeReads(input, read_base_index, hap_base_index, host_input[0]->numRead, host_input[1]->numRead, host_input[2]->numRead, cur_numRead, cur_numHap, maxCols, violate);
    //fprintf(stderr, "numRead0 = %d, numRead1 = %d, numRead2 = %d, total read = %d\n", host_input[0]->numRead, host_input[1]->numRead, host_input[2]->numRead, cur_numRead);
    if(violate){
        //fprintf(stderr, "Either one SLR has 0 reads or this segment is too small to be run on FPGA, switch back to AVX\n");
        return 0;
    }
    sortReads(input, read_base_index, cur_numRead, cur_numHap, maxCols, slr_pu_num, host_input);
 
    return cells;
}

JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_computeLikelihoodsNative
(JNIEnv* env, jobject obj, jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray){
    struct PairhmmContext* curCtx = init_pairhmm_context();
    struct timespec time_start, time_end, diff_overall;
    struct timespec time1, time2, diff;
    clock_gettime(CLOCK_REALTIME, &time_start);

    //==================================================================
    // get Java data
    //printf("in the native, start to get the input data\n");
    int numRead = env->GetArrayLength(readDataArray);
    int numHap = env->GetArrayLength(haplotypeDataArray);
    bool useFPGA = !curCtx->g_use_double;
    JavaData javaData;
    javaData.getData(env, readDataArray, haplotypeDataArray, curCtx->client_input, numRead, numHap, useFPGA);
    curCtx->javaResults = javaData.getOutputArray(env, likelihoodArray);
    if(!useFPGA){
        compute_violate(curCtx->client_input, curCtx->javaResults, curCtx->g_use_double);
        javaData.releaseData(env);
        clock_gettime(CLOCK_REALTIME, &time_end);
        diff_overall = pmm_diff_time(time_start, time_end);
        curCtx->total_time += diff_overall.tv_sec * 1e9 + diff_overall.tv_nsec;
        return;
    }
    curCtx->host_raw_output.clear();
    curCtx->host_raw_output.resize(numRead * numHap);
    int max_rsdata_num = MAX_RSDATA_NUM;
    int max_hapdata_num = MAX_HAPDATA_NUM;
    curCtx->pmm_count++;
    int numReadSeg = numRead / max_rsdata_num - (numRead % max_rsdata_num == 0) + 1;
    int numHapSeg = numHap / max_hapdata_num - (numHap % max_hapdata_num == 0) + 1;
    bool violate = false;
    long long base = 0;
    curCtx->read_base_index = 0;
    curCtx->hap_base_index = 0;
    for(int i = 0; i < numReadSeg; i++){
        curCtx->cur_numRead = i + 1 < numReadSeg ? max_rsdata_num : numRead - i * max_rsdata_num;
        curCtx->hap_base_index = 0;
        for(int j = 0; j < numHapSeg; j++){
            base = curCtx->read_base_index * numHap + curCtx->hap_base_index; 
            curCtx->cur_numHap = j + 1 < numHapSeg ? max_hapdata_num: numHap - j * max_hapdata_num;
            float cells = update_host_inputs(curCtx->client_input, curCtx->read_base_index, curCtx->hap_base_index, curCtx->cur_numRead, curCtx->cur_numHap, violate, curCtx->host_input);
            if(violate){
                computePairhmmAVXSegment(curCtx->client_input, curCtx->read_base_index, curCtx->hap_base_index, curCtx->cur_numRead, curCtx->cur_numHap, curCtx->host_raw_output);
            }
            else{
                clock_gettime(CLOCK_REALTIME, &time1);
                PairhmmClient client;
                curCtx->rejected = false;
                //fprintf(stderr, "InputDataPackOpt has %d bytes, %d 512-bit %f MB\n", sizeof(InputDataPackOpt), sizeof(InputDataPackOpt) / 64, (double)sizeof(InputDataPackOpt) / 1024 / 1024);
                client.setInput(0, (char*)(&(curCtx->host_input[0]->dataPack)), inputBank[0], 1, sizeof(InputDataPackOpt), sizeof(char), BLAZE_INPUT);
                client.setInput(1, &(curCtx->host_input[0]->numRead), 1, 1, sizeof(int), BLAZE_INPUT);
                client.setInput(2, &(curCtx->host_input[0]->numHap), 1, 1, sizeof(int), BLAZE_INPUT);
                client.setInput(3, (char*)(&(curCtx->host_input[1]->dataPack)), inputBank[1], 1, sizeof(InputDataPackOpt), sizeof(char), BLAZE_INPUT);
                client.setInput(4, &(curCtx->host_input[1]->numRead), 1, 1, sizeof(int), BLAZE_INPUT);
                client.setInput(5, &(curCtx->host_input[1]->numHap), 1, 1, sizeof(int), BLAZE_INPUT);
#if DIE_NUM == 3
                client.setInput(6, (char*)(&(curCtx->host_input[2]->dataPack)), inputBank[2], 1, sizeof(InputDataPackOpt), sizeof(char), BLAZE_INPUT);
                client.setInput(7, &(curCtx->host_input[2]->numRead), 1, 1, sizeof(int), BLAZE_INPUT);
                client.setInput(8, &(curCtx->host_input[2]->numHap), 1, 1, sizeof(int), BLAZE_INPUT);
#endif
                client.setInput(3 * DIE_NUM, &(cells), 1, 1, sizeof(float), BLAZE_INPUT);
                clock_gettime(CLOCK_REALTIME, &time2);
                diff = pmm_diff_time(time1, time2);
                curCtx->fpga_prepare_data_time += diff.tv_sec * 1e9 + diff.tv_nsec;


                clock_gettime(CLOCK_REALTIME, &time1);
                client.start();
                if(!curCtx->rejected){
                    float* raw_output0 = (float*)client.getOutputPtr(0);
                    float* raw_output1 = (float*)client.getOutputPtr(1);
#if DIE_NUM == 3
                    float* raw_output2 = (float*)client.getOutputPtr(2);
#endif
                    for(int k = 0; k < curCtx->cur_numHap * curCtx->cur_numRead; k++){
                        float cur_float_result = 0;
                        if(k < curCtx->host_input[0]->numRead * curCtx->cur_numHap)
                            cur_float_result = raw_output0[k];
                        else if(k - curCtx->host_input[0]->numRead * curCtx->cur_numHap < curCtx->host_input[1]->numRead * curCtx->cur_numHap) 
                            cur_float_result = raw_output1[k - curCtx->host_input[0]->numRead * curCtx->cur_numHap];
#if DIE_NUM == 3
                        else
                            cur_float_result = raw_output2[k - curCtx->host_input[0]->numRead * curCtx->cur_numHap - curCtx->host_input[1]->numRead * curCtx->cur_numHap];
#endif
                        curCtx->host_raw_output[base + (k / curCtx->cur_numHap) * numHap + (k % curCtx->cur_numHap)] = cur_float_result;
                    }
                }
                clock_gettime(CLOCK_REALTIME, &time2);
                diff = pmm_diff_time(time1, time2);
                curCtx->fpga_time += diff.tv_sec * 1e9 + diff.tv_nsec;
            }
            violate = false;

            curCtx->hap_base_index += curCtx->cur_numHap;
        }
        curCtx->read_base_index += curCtx->cur_numRead;
    }
    
    clock_gettime(CLOCK_REALTIME, &time1);
    for(int k = 0; k < numRead * numHap; k++){
        if(curCtx->host_raw_output[k] < MIN_ACCEPTED){
            testcase testCase;
            int ii = k / numHap;
            int jj = k % numHap;
            testCase.rslen = curCtx->client_input->reads[ii].bases.size();
            testCase.rs = curCtx->client_input->reads[ii].bases.c_str();
            testCase.i = curCtx->client_input->reads[ii]._i.c_str();
            testCase.d = curCtx->client_input->reads[ii]._d.c_str();
            testCase.q = curCtx->client_input->reads[ii]._q.c_str();
            testCase.c = curCtx->client_input->reads[ii]._c.c_str();
            testCase.haplen = curCtx->client_input->haps[jj].bases.size();
            testCase.hap = curCtx->client_input->haps[jj].bases.c_str();
            double result_double = compute_fp_avxd(&testCase);
            curCtx->javaResults[k] = log10(result_double) - g_ctxd.LOG10_INITIAL_CONSTANT;  
        }
        else{
            curCtx->javaResults[k] = log10(curCtx->host_raw_output[k]) - g_ctxf.LOG10_INITIAL_CONSTANT;
        }
    }
    clock_gettime(CLOCK_REALTIME, &time2);
    diff = pmm_diff_time(time1, time2);
    curCtx->overflow_avx_time += diff.tv_sec * 1e9 + diff.tv_nsec;

    //==================================================================
    // release Java data
    javaData.releaseData(env);
    
    clock_gettime(CLOCK_REALTIME, &time_end);
    diff_overall = pmm_diff_time(time_start, time_end);
    curCtx->total_time += diff_overall.tv_sec * 1e9 + diff_overall.tv_nsec;
    
}



JNIEXPORT void JNICALL Java_com_falconcomputing_genomics_haplotypecaller_FalconPairhmm_doneNative
(JNIEnv* env, jobject obj)
{
    struct PairhmmContext* curCtx = init_pairhmm_context();
    DLOG(INFO) << "pmm_call_count = " << curCtx->pmm_count << ", total native time = " << curCtx->total_time * 1e-9 << " secs";
    DLOG(INFO) << "fpga time = " << curCtx->fpga_time * 1e-9 << " secs";
    DLOG(INFO) << "fpga prepare time = " << curCtx->fpga_prepare_data_time * 1e-9 << " secs";
    DLOG(INFO) << "directx avx time = " << curCtx->direct_avx_time * 1e-9 << " secs";
    DLOG(INFO) << "overflow recompute avx time = " << curCtx->overflow_avx_time * 1e-9 << " secs";
    DLOG(INFO) << "reject avx time = " <<  curCtx->reject_time * 1e-9 << " secs";
    delete curCtx->client_input;
    _mm_free(curCtx->host_input[0]);
    _mm_free(curCtx->host_input[1]);
    _mm_free(curCtx->host_input[2]);
    curCtx->host_raw_output.clear();
    delete curCtx;
    std::thread::id this_id = std::this_thread::get_id();
    pmm_ctx_map_mutex.lock();
    pmm_map.erase(this_id);   
    pmm_ctx_map_mutex.unlock();
}
