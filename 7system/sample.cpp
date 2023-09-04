#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define CHROM_LENGTH 500 // ��`�q�̒����A���i���Ɠ���
#define GENERATION_MAX		100	// �����㐔
#define POP_SIZE			100	// �̏W�c���̌̂̐�
#define ELITE 1

//�{�Ԃ̕]���֐��ł����L�Ɠ����l���g���܂�
#define CAL_LIMIT_PER_DAY 2400 //������̐H���̐ێ�J�����[�̏���A����̋�؂�����߂�̂Ɏg�p�A���Ō��߂Ă���̂ō����͖���
#define IDEAL_CAL_PER_DAY 2164 //���������̗��z�I�Ȑڎ�J�����[

// 0�ȏ�1�ȉ��̎�������
#define RAND_01 ((double)rand() / RAND_MAX)
// �W���̗v�f�ƂȂ�ő吔�̕����l
#define N				64		
//???
#define DBL_MAX         1.7976931348623158e+308
// �g�[�i�����g�T�C�Y
#define TOURNAMENT_SIZE	30		

// �K���x��ϊ������l
double trFit[POP_SIZE];

typedef struct individual {
    double fitness;			// �K���x
    int chrom[CHROM_LENGTH];		// ���F��
} Individual;

typedef struct {
    int id;
    //char name[1024]
    double calorie;
    double protein;
    double fat;
    double carbohydrate;
}Product;

void load_products_data(char* fileName, Product* products) {
    int i, columnNum;
    char line[1024];
    FILE* fp;
    char* pos1;
    char* pos2;

    int dataNum;      // �f�[�^��
    int colNum;     // ���ڐ�

    // �����ϐ��̐��ƃf�[�^���̎擾
    if ((fp = fopen(fileName, "r")) == NULL) {
        printf("%s���J���܂���D\n", fileName);
        exit(1);
    }
    colNum = -1;
    dataNum = 0;
    while (fgets(line, 1024, fp)) {
        if (strcmp(line, "\n")) {
            columnNum = 1;
            pos1 = line;
            do {
                pos2 = strchr(pos1, ',');
                if (pos2) {
                    if (pos2 == pos1) {
                        printf("��̃f�[�^���܂܂�Ă��܂��D");
                        exit(1);
                    }
                    columnNum++;
                    pos1 = pos2 + 1;
                }
            } while (pos2);
            if (*pos1 == '\n') {
                printf("��̃f�[�^���܂܂�Ă��܂��D");
                exit(1);
            }
            if (colNum == -1) {
                colNum = columnNum - 1;
            }
            else if (colNum != columnNum - 1) {
                printf("�񐔂̈قȂ郌�R�[�h������܂��D");
                exit(1);
            }
            dataNum++;
        }
    }
    fclose(fp);

    // �f�[�^��Ǎ���
    if ((fp = fopen(fileName, "r")) == NULL) {
        printf("%s���J���܂���D\n", fileName);
        exit(1);
    }
    fgets(line, 1024, fp);

    // �w�b�_�[��ǂݔ�΂�
    pos1 = line;
    pos2 = strchr(pos1, '\n');
    pos1 = pos2 + 1;

    for (i = 0; i < dataNum - 1; i++) {
        fgets(line, 1024, fp);
        pos1 = line;
        //printf("%s", line);
        // id
        pos2 = strchr(pos1, ',');
        *pos2 = '\0';
        //printf("%s\n", pos1);
        products[i].id = atoi(pos1);
        pos1 = pos2 + 1;

        //calorie
        pos2 = strchr(pos1, ',');
        *pos2 = '\0';
        products[i].calorie = atof(pos1);
        pos1 = pos2 + 1;

        //protein
        pos2 = strchr(pos1, ',');
        *pos2 = '\0';
        products[i].protein = atof(pos1);
        pos1 = pos2 + 1;

        //fat
        pos2 = strchr(pos1, ',');
        *pos2 = '\0';
        products[i].fat = atof(pos1);
        pos1 = pos2 + 1;

        //carbohydrate
        pos2 = strchr(pos1, '\n');
        *pos2 = '\0';
        products[i].carbohydrate = atof(pos1);
        pos1 = pos2 + 1;
    }
    fclose(fp);
    /*
    for(i = 0; i < dataNum -1; i++){
       printf("id: %d\n\tcalotie: %f\n\tp: %f\n\tf: %f\n\tc: %f\n",
       products[i].id,
       products[i].calorie,
       products[i].protein,
       products[i].fat,
       products[i].carbohydrate);
    }*/

}

//��`�q�̓p�X�\���ŁA��`�q������
void init_chrom(int* chrom) {
    int i, j, k;
    int tmp;
    for (int i = 0; i < CHROM_LENGTH; i++) {
        chrom[i] = i;
    }
    //�V���b�t��
    for (i = 0; i < CHROM_LENGTH; i++) {
        j = rand() % CHROM_LENGTH;
        k = rand() % CHROM_LENGTH;
        tmp = chrom[j];
        chrom[j] = chrom[k];
        chrom[k] = tmp;
    }

    /*
    printf("init\n");
    printf("DEBUG:\n[");
    for(int i = 0; i < CHROM_LENGTH; i++){
        printf("%d, ", chrom[i]);
    }
    printf("]\n");
    */

}

int exchange(int* p1, int* p2, int len, int value) {
    int i;
    for (i = 0; i < len; i++) {
        if (value == p1[i]) {
            //printf("%d -> ",value);
            return exchange(p1, p2, len, p2[i]);
        }
    }
    //printf("%d\n",value);
    return value;
}

//�����ʑ�����
void crossoverPMX(int* p1, int* p2, int* c1, int* c2) {
    int i;
    int index1, index2, body_length;
    //�����_���I��
    index1 = rand() % CHROM_LENGTH;
    body_length = rand() % (CHROM_LENGTH - index1);//�K�v�ȗ��R
    index2 = index1 + body_length;

    //c1 ��1��(0 - index1-1 )
    for (i = 0; i < index1; i++) {
        c1[i] = exchange(&p2[index1], &p1[index1], body_length, p1[i]);
    }
    //c1 ��2��(index1 - index2-1)
    for (i = index1; i < index2; i++) {
        c1[i] = p2[i];
    }
    //c1 ��3�� (index2 - )
    for (i = index2; i < CHROM_LENGTH; i++) {
        c1[i] = exchange(&p2[index1], &p1[index1], body_length, p1[i]);
    }

    //c2 ��1��(0 - index1-1 )
    for (i = 0; i < index1; i++) {
        c2[i] = exchange(&p1[index1], &p2[index1], body_length, p2[i]);
    }
    //c2 ��2��(index1 - index2-1)
    for (i = index1; i < index2; i++) {
        c2[i] = p1[i];
    }
    //c2 ��3�� (index2 - )
    for (i = index2; i < CHROM_LENGTH; i++) {
        c2[i] = exchange(&p1[index1], &p2[index1], body_length, p2[i]);
    }

}

//��������
void crossoverOX(int* p1, int* p2, int* c1, int* c2) {
    int i;
    int index1, index2, body_length;

    //c1 ��1��(0 - index1-1 )
    for (i = 0; i < index1; i++) {
        c1[i] = exchange(&p2[index1], &p1[index1], body_length, p1[i]);
    }
    //c1 ��2���@(index1 - index2-1)
    for (i = index1; i < index2; i++) {
        c1[i] = p2[i];
    }
}


double calc_total_calorie(Product* products) {
    double total_calorie;
    total_calorie = 0;
    for (int i = 0; i < CHROM_LENGTH; i++) {
        total_calorie += products[i].calorie;
    }
    return total_calorie;
}

void mutate(int* chrom) {}  //�ق���������ǉ�

// �����L���O�I��
int rankingSelect1()
{
    int num, denom, r;

    denom = POP_SIZE * (POP_SIZE + 1) / 2;
    r = ((rand() << 16) + (rand() << 1) + (rand() % 2)) % denom + 1;
    for (num = POP_SIZE; 0 < num; num--) {
        if (r <= num) {
            break;
        }
        r -= num;
    }
    return POP_SIZE - num;
}

// �m���Ɋ�Â������L���O�I���Őe�̂�1�I������
// �߂�l: �I�������e�̂̓Y����
int rankingSelect2()
{
    int rank, denom;
    double prob, r;

    denom = POP_SIZE * (POP_SIZE + 1) / 2;
    r = RAND_01;
    for (rank = 1; rank < POP_SIZE; rank++) {
        prob = (double)(POP_SIZE - rank + 1) / denom;
        if (r <= prob) {
            break;
        }
        r -= prob;
    }
    return rank - 1;
}

// ���[���b�g�I���Őe�̂�1�I������
// �߂�l: �I�������e�̂̓Y����
int rouletteSelect(double* trFit, double denom)
{
    int rank;
    double prob, r;


    r = RAND_01;
    for (rank = 1; rank < POP_SIZE; rank++) {
        prob = trFit[rank - 1] / denom;
        if (r <= prob) {
            break;
        }
        r -= prob;
    }
    return rank - 1;
}

// �g�[�i�����g�I���Őe�̂�1�I������
// �߂�l: �I�������e�̂̓Y����
int tournamentSelect(Individual* population)
{
    int i, ret, num, r;
    double bestFit;
    int tmp[N];

    for (i = 0; i < N; i++) {
        tmp[i] = 0;
    }
    ret = -1;
    bestFit = DBL_MAX;
    num = 0;
    while (1) {
        r = rand() % N;
        if (tmp[r] == 0) {
            tmp[r] = 1;
            if (population[r].fitness < bestFit) {
                ret = r;
                bestFit = population[r].fitness;
            }
            if (++num == TOURNAMENT_SIZE) {
                break;
            }
        }
    }
    return ret;
}

void new_generation(Individual* population) {
    int i;
    int parent_index1, parent_index2;
    // ���[���b�g�I���̊m�������߂�Ƃ��̕���
    double denom = 0.0;		
	

    Individual* next_pop;
    next_pop = (Individual*)malloc(POP_SIZE * sizeof(Individual));

    // ���[���b�g�I���̂��߂̏���
    /*
    for(i = 0; i < POP_SIZE; i++) {
        trFit[i] = (population[POP_SIZE - 1].fitness - population[i].fitness) / (population[POP_SIZE - 1].fitness - population[0].fitness);
        denom += trFit[i];
    }
    */

    parent_index1 = 0;
    parent_index2 = 1;
    for (i = ELITE; i < POP_SIZE - 1; i += 2) {
        //�e�I��
        // ���ʂɊ�Â������L���O�I��
        parent_index1 = rankingSelect1();
        parent_index2 = rankingSelect1();
        // �m���Ɋ�Â������L���O�I��
        //parent_index1 = rankingSelect2();
        //parent_index2 = rankingSelect2();
        // ���[���b�g�I��
        //parent_index1 = rouletteSelect(trFit,denom);
        //parent_index2 = rouletteSelect(trFit,denom);
        // �g�[�i�����g�I��
        //parent_index1 = tournamentSelect(population);
        //parent_index2 = tournamentSelect(population);
        
        //�����ʑ�����
        crossoverPMX(population[parent_index1].chrom, population[parent_index2].chrom, next_pop[i].chrom, next_pop[i + 1].chrom);
        //��������
        //crossoverOX(population[parent_index1].chrom, population[parent_index2].chrom, next_pop[i].chrom, next_pop[i + 1].chrom);
    }
    //��z�ɂ���č��̂̐�����̎��A�Ō�̗]�����g�̓����_���ɐ��������̂�����
    if (i != POP_SIZE) {
        init_chrom(next_pop[i].chrom);
    }
    //�ˑR�ψ�
    for (i = ELITE; i < POP_SIZE; i++) {
        mutate(next_pop[i].chrom);
    }
    for (i = ELITE; i < POP_SIZE; i++) {
        population[i] = next_pop[i];
    }
    free(next_pop);

}

// �z��population[lb]�`population[ub]��fitness�̏����ɐ���i�N�C�b�N�\�[�g�j
// lb: ���񂷂�͈͂̉����̓Y��
// ub: ���񂷂�͈͂̏���̓Y��
void sort(Individual* population, int lb, int ub)
{
    int i, j, k;
    double pivot;
    Individual tmp;

    if (lb < ub) {
        k = (lb + ub) / 2;
        pivot = population[k].fitness;
        i = lb;
        j = ub;
        do {
            while (population[i].fitness < pivot)
                i++;
            while (population[j].fitness > pivot)
                j--;
            if (i <= j) {
                tmp = population[i];
                population[i] = population[j];
                population[j] = tmp;
                i++;
                j--;
            }
        } while (i <= j);
        sort(population, lb, j);
        sort(population, i, ub);
    }
}

//�]���֐�
void evaluation(Product* products, Individual* population) {
    int i, j;
    int term; //����
    int product_index;
    //�^���p�N���A�����A�Y�������A�G�l���M�[�̊e���ɂ����鍇�v��
    double day_protein, day_fat, day_carb, day_calorie;
    //�^���p�N���A�����A�Y�������A�G�l���M�[�A�e���̈�E�x
    double protein_deviance, fat_deviance, carb_deviance, calorie_deviance, day_deviance;
    double ave, sd, sqrSum;
    for (i = 0; i < POP_SIZE; i++) {
        j = 0;
        term = 0;
        day_protein = day_fat = day_carb = day_calorie = 0.0;
        ave = 0.0;
        sqrSum = 0.0;
        while (j < CHROM_LENGTH) {
            product_index = population[i].chrom[j];
            // printf("id:%d\n", products[product_index].id);
            // printf("protein:%f\n", products[product_index].protein);

            if (day_calorie < CAL_LIMIT_PER_DAY) {
                day_protein += products[product_index].protein;
                day_fat += products[product_index].fat;
                day_carb += products[product_index].carbohydrate;
                day_calorie += products[product_index].calorie;
                j++;
            }
            else {
                protein_deviance = fabs(0.2 - (double)(day_protein * 4) / day_calorie) * (1 / 0.2);
                fat_deviance = fabs(0.2 - (double)(day_fat * 9) / day_calorie) * (1 / 0.2);
                carb_deviance = fabs(0.6 - (double)(day_carb * 9) / day_calorie) * (1 / 0.6);
                calorie_deviance = fabs(IDEAL_CAL_PER_DAY - day_calorie) * (1 / IDEAL_CAL_PER_DAY);
                day_deviance = exp(protein_deviance) + exp(fat_deviance) + exp(carb_deviance) + exp(calorie_deviance) - 4;

                //�Ō�ɕ��ϒl�ƕW���΍������߂邽�߂ɑ����Ă���
                ave += day_deviance;
                sqrSum += day_deviance * day_deviance;
                term++;
                // printf("ave:%f\n", ave);
                // printf("sqrSum:%f\n", sqrSum);

                day_calorie = 0.0;
                day_protein = 0.0;
                day_fat = 0.0;
                day_carb = 0.0;
            }
        }
        //�ŏI���̈�E�x���v�Z
        if (day_calorie != 0.0) {
            protein_deviance = fabs(0.2 - (double)(day_protein * 4) / day_calorie) * (1 / 0.2);
            fat_deviance = fabs(0.2 - (double)(day_fat * 9) / day_calorie) * (1 / 0.2);
            carb_deviance = fabs(0.6 - (double)(day_carb * 9) / day_calorie) * (1 / 0.6);
            calorie_deviance = fabs(IDEAL_CAL_PER_DAY - day_calorie) * (1 / IDEAL_CAL_PER_DAY);
            day_deviance = exp(protein_deviance) + exp(fat_deviance) + exp(carb_deviance) + exp(calorie_deviance) - 4;

            ave += day_deviance;
            sqrSum += day_deviance * day_deviance;
            term++;
        }

        ave /= term;                           //��E�x�̕��ϒl
        sd = sqrt(sqrSum / term - ave * ave);  //��E�x�̕W���΍�
        //printf("term:%d\n", term);
        //printf("innerdef_ave:%f\n", ave);
        // printf("innerdef_sqrSum:%f\n", sqrSum);
        //printf("innerdef_sd:%f\n", sd);
        population[i].fitness = ave + sd;
        //printf("innerdef_fitness:%f\n", population[i].fitness);
    }
    sort(population, 0, POP_SIZE - 1);
}

void save_result(char* fileName, int* chrom) {
    FILE* fp;
    int i;
    // ������
    if ((fp = fopen(fileName, "w")) == NULL) {
        printf("�t�@�C��%s���J���܂���D\n", fileName);
        exit(-1);
    }
    for (i = 0; i < CHROM_LENGTH; i++) {
        fprintf(fp, "%d,", chrom[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

int main() {  
    int i;
    int gen;
    Product products[CHROM_LENGTH];
    char fname[64];

    Individual* population;
    population = (Individual*)malloc(POP_SIZE * sizeof(Individual));

    // �����̃^�l�̐ݒ�
    srand((unsigned int)time(NULL));

    //�f�[�^��products�ɓǂݍ���
    strcpy(fname, "dev_selected_500.csv");
    load_products_data(fname, products);

    // �����W�c����
    for (i = 0; i < POP_SIZE; i++) {
        init_chrom(population[i].chrom);
    }

    // �i��
    for (gen = 1; gen <= GENERATION_MAX; gen++) {
        // �����㐶��
        new_generation(population);
        // �]��
        evaluation(products, population);
        // �r���o�ߕ\��
        fprintf((__acrt_iob_func(2)), "��%4d����\t�ŗD�ǌ̂̓K���x��%.16lf\n", gen, population[0].fitness);
    }
    strcpy(fname, "result.csv");
    save_result(fname, population[0].chrom);
    free(population);
    return 0;
}