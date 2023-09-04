#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define CHROM_LENGTH 500 // 遺伝子の長さ、商品数と同じ
#define GENERATION_MAX		100	// 世代交代数
#define POP_SIZE			100	// 個体集団内の個体の数
#define ELITE 1

//本番の評価関数でも下記と同じ値を使います
#define CAL_LIMIT_PER_DAY 2400 //一日分の食事の摂取カロリーの上限、一日の区切りを決めるのに使用、勘で決めているので根拠は無い
#define IDEAL_CAL_PER_DAY 2164 //一日あたりの理想的な接種カロリー

// 0以上1以下の実数乱数
#define RAND_01 ((double)rand() / RAND_MAX)
// 集合の要素となる最大数の平方値
#define N				64		
//???
#define DBL_MAX         1.7976931348623158e+308
// トーナメントサイズ
#define TOURNAMENT_SIZE	30		

// 適応度を変換した値
double trFit[POP_SIZE];

typedef struct individual {
    double fitness;			// 適応度
    int chrom[CHROM_LENGTH];		// 染色体
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

    int dataNum;      // データ数
    int colNum;     // 項目数

    // 説明変数の数とデータ数の取得
    if ((fp = fopen(fileName, "r")) == NULL) {
        printf("%sが開けません．\n", fileName);
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
                        printf("空のデータが含まれています．");
                        exit(1);
                    }
                    columnNum++;
                    pos1 = pos2 + 1;
                }
            } while (pos2);
            if (*pos1 == '\n') {
                printf("空のデータが含まれています．");
                exit(1);
            }
            if (colNum == -1) {
                colNum = columnNum - 1;
            }
            else if (colNum != columnNum - 1) {
                printf("列数の異なるレコードがあります．");
                exit(1);
            }
            dataNum++;
        }
    }
    fclose(fp);

    // データを読込む
    if ((fp = fopen(fileName, "r")) == NULL) {
        printf("%sが開けません．\n", fileName);
        exit(1);
    }
    fgets(line, 1024, fp);

    // ヘッダーを読み飛ばす
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

//遺伝子はパス表現で、遺伝子初期化
void init_chrom(int* chrom) {
    int i, j, k;
    int tmp;
    for (int i = 0; i < CHROM_LENGTH; i++) {
        chrom[i] = i;
    }
    //シャッフル
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

//部分写像交叉
void crossoverPMX(int* p1, int* p2, int* c1, int* c2) {
    int i;
    int index1, index2, body_length;
    //交叉点を二つ選択
    index1 = rand() % CHROM_LENGTH;
    body_length = rand() % (CHROM_LENGTH - index1);//必要な理由
    index2 = index1 + body_length;

    //c1 第1部(0 - index1-1 )
    for (i = 0; i < index1; i++) {
        c1[i] = exchange(&p2[index1], &p1[index1], body_length, p1[i]);
    }
    //c1 第2部(index1 - index2-1)
    for (i = index1; i < index2; i++) {
        c1[i] = p2[i];
    }
    //c1 第3部 (index2 - )
    for (i = index2; i < CHROM_LENGTH; i++) {
        c1[i] = exchange(&p2[index1], &p1[index1], body_length, p1[i]);
    }

    //c2 第1部(0 - index1-1 )
    for (i = 0; i < index1; i++) {
        c2[i] = exchange(&p1[index1], &p2[index1], body_length, p2[i]);
    }
    //c2 第2部(index1 - index2-1)
    for (i = index1; i < index2; i++) {
        c2[i] = p1[i];
    }
    //c2 第3部 (index2 - )
    for (i = index2; i < CHROM_LENGTH; i++) {
        c2[i] = exchange(&p1[index1], &p2[index1], body_length, p2[i]);
    }

}

//順序交叉
void crossoverOX(int* p1, int* p2, int* c1, int* c2) {
    int i;
    int index1, index2, body_length;

    //c1 第1部(0 - index1-1 )
    for (i = 0; i < index1; i++) {
        c1[i] = exchange(&p2[index1], &p1[index1], body_length, p1[i]);
    }
    //c1 第2部　(index1 - index2-1)
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

void mutate(int* chrom) {}  //ほしかったら追加

// ランキング選択
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

// 確率に基づくランキング選択で親個体を1つ選択する
// 戻り値: 選択した親個体の添え字
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

// ルーレット選択で親個体を1つ選択する
// 戻り値: 選択した親個体の添え字
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

// トーナメント選択で親個体を1つ選択する
// 戻り値: 選択した親個体の添え字
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
    // ルーレット選択の確率を求めるときの分母
    double denom = 0.0;		
	

    Individual* next_pop;
    next_pop = (Individual*)malloc(POP_SIZE * sizeof(Individual));

    // ルーレット選択のための処理
    /*
    for(i = 0; i < POP_SIZE; i++) {
        trFit[i] = (population[POP_SIZE - 1].fitness - population[i].fitness) / (population[POP_SIZE - 1].fitness - population[0].fitness);
        denom += trFit[i];
    }
    */

    parent_index1 = 0;
    parent_index2 = 1;
    for (i = ELITE; i < POP_SIZE - 1; i += 2) {
        //親選択
        // 順位に基づくランキング選択
        parent_index1 = rankingSelect1();
        parent_index2 = rankingSelect1();
        // 確率に基づくランキング選択
        //parent_index1 = rankingSelect2();
        //parent_index2 = rankingSelect2();
        // ルーレット選択
        //parent_index1 = rouletteSelect(trFit,denom);
        //parent_index2 = rouletteSelect(trFit,denom);
        // トーナメント選択
        //parent_index1 = tournamentSelect(population);
        //parent_index2 = tournamentSelect(population);
        
        //部分写像交叉
        crossoverPMX(population[parent_index1].chrom, population[parent_index2].chrom, next_pop[i].chrom, next_pop[i + 1].chrom);
        //順序交叉
        //crossoverOX(population[parent_index1].chrom, population[parent_index2].chrom, next_pop[i].chrom, next_pop[i + 1].chrom);
    }
    //交配によって作る個体の数が奇数の時、最後の余った枠はランダムに生成した個体を入れる
    if (i != POP_SIZE) {
        init_chrom(next_pop[i].chrom);
    }
    //突然変異
    for (i = ELITE; i < POP_SIZE; i++) {
        mutate(next_pop[i].chrom);
    }
    for (i = ELITE; i < POP_SIZE; i++) {
        population[i] = next_pop[i];
    }
    free(next_pop);

}

// 配列population[lb]～population[ub]をfitnessの昇順に整列（クイックソート）
// lb: 整列する範囲の下限の添字
// ub: 整列する範囲の上限の添字
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

//評価関数
void evaluation(Product* products, Individual* population) {
    int i, j;
    int term; //日数
    int product_index;
    //タンパク質、脂質、炭水化物、エネルギーの各日における合計量
    double day_protein, day_fat, day_carb, day_calorie;
    //タンパク質、脂質、炭水化物、エネルギー、各日の逸脱度
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

                //最後に平均値と標準偏差を求めるために足しておく
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
        //最終日の逸脱度を計算
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

        ave /= term;                           //逸脱度の平均値
        sd = sqrt(sqrSum / term - ave * ave);  //逸脱度の標準偏差
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
    // 書込み
    if ((fp = fopen(fileName, "w")) == NULL) {
        printf("ファイル%sが開けません．\n", fileName);
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

    // 乱数のタネの設定
    srand((unsigned int)time(NULL));

    //データをproductsに読み込む
    strcpy(fname, "dev_selected_500.csv");
    load_products_data(fname, products);

    // 初期集団生成
    for (i = 0; i < POP_SIZE; i++) {
        init_chrom(population[i].chrom);
    }

    // 進化
    for (gen = 1; gen <= GENERATION_MAX; gen++) {
        // 次世代生成
        new_generation(population);
        // 評価
        evaluation(products, population);
        // 途中経過表示
        fprintf((__acrt_iob_func(2)), "第%4d世代\t最優良個体の適応度は%.16lf\n", gen, population[0].fitness);
    }
    strcpy(fname, "result.csv");
    save_result(fname, population[0].chrom);
    free(population);
    return 0;
}