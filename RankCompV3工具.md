## RankCompV3.jl

**严婧**

**Yan Jing**

**yanjer123@qq.com**



​	基于相对表达秩序关系REO设计的识别差异表达基因的软件包。工具基于julia语言开发，已配备软件包可供直接使用，具体介绍如下。

### 1、RankCompV3适用范围

| 分类                   | 适用范围                                                     |
| ---------------------- | ------------------------------------------------------------ |
| **物种**               | 不限定物种                                                   |
| **管家基因支持的物种** | human/mouse                                                  |
| **数据类型**           | 芯片数据、高通量测序数据（RNA-seq）、单细胞测序数据（scRNA-seq）、CEL-seq测序数据、甲基化数据和蛋白质组数据 |
| **数据格式**           | Count/log2Count/FPKM/RPKM/TPM/表达值                         |
| **基因类型**           | REFSEQ/SYMBOL/ENTREZID/ENSEMBL/UNIGENE/GENENAME              |

### 2、输入文件说明

| 文件                         | 格式                                                         |
| ---------------------------- | ------------------------------------------------------------ |
| **Expression profile file ** | Each row represents the gene, each column represents the sample and the expression matrix of the gene in the first column. |
| **Metadata file **           | First column sample name, second column group information.   |

### 3、工具依赖package

​	**其中，julia相关软件包已加入环境，无需下载。R语言相关软件包会进行自动下载，您也可手动进行下载。**

#### 1）Julia:

​	Distributed, SharedArrays, Base.Threads, MultipleTesting, DataFrames, DelimitedFiles, CSV, Statistics, Rcall, ArgParse

#### 2）R语言:

​	ggplot2, pheatmap

### 3、主函数reoa()参数说明

| 参数名              | 参数类型                     | default           | 参数描述                                                     |
| ------------------- | ---------------------------- | ----------------- | ------------------------------------------------------------ |
| fn_expr             | AbstractString               |                   | Gene expression profile file path.（required）               |
| fn_metadata         | AbstractStringAbstractString |                   | Grouping information file path.（required）                  |
| expr_threshold      | NumberNumber                 | 3                 | Gene expression threshold.                                   |
| pval_reo            | AbstractFloatAbstractFloat   | 0.01              | Stable threshold for p-value.                                |
| pval_sign_reo       | AbstractFloatAbstractFloat   | 1                 | Significant  reversal threshold for p-value.                 |
| padj_sign_reo       | AbstractFloatAbstractFloat   | 0.05              | Significant reversal threshold for FDR  value.               |
| hk_file             | AbstractString               | HK_genes_info.tsv | Housekeeper gene  file path.                                 |
| hk_name             | AbstractString               | ENSEMBL           | Column name of the  column where the housekeeping gene is located. |
| ref_gene_num        | Int                          | 3000              | The upper limit of  the number of housekeeping genes, if it is greater than this value,  ref_gene_num housekeeping genes are randomly selected from it. |
| no_use_housekeeping | Int                          | 0                 | Do not use  housekeeping gene to set 1, use housekeeping gene to set 0. |
| species             | AbstractString               | human             | Species information of the data,  currently supports hunan and mouse. |
| cell_drop_rate      | Int                          | 0                 | At least how many genes were detected in each sample. (value non-zero). |
| gene_drop_rate      | Int                          | 0                 | At least how many cells each gene was  detected in (value non-zero). |
| work_dir            | AbstractString               | ./                | Working Directory.                                           |

### 4、工具使用方式

#### 1）julia中RankCompV3包

[Yanjj1/RankCompV3: Differential expression gene recognition algorithm based on REOs. (github.com)](https://github.com/Yanjj1/RankCompV3)

##### （1）直接使用

```julia
#在julia中载入Pkg
using Pkg
#首次使用需要安装RankCompV3包进行使用
Pkg.add(url="https://github.com/Yanjj1/RankCompV3.git/master")
using RankCompV3
#runing
RankCompV3.reoa(fn_expr = "expr.txt",
        fn_metadata = "metadata.txt";
    expr_threshold = 3,
          pval_reo = 0.01,
     pval_sign_reo = 1.00,
     padj_sign_reo = 0.05,
           hk_file = "HK_genes_info.tsv",
           hk_name = "ENSEMBL",
      ref_gene_num = 3000,
    no_use_housekeeping = 0,
           species = "human",
    cell_drop_rate = 0,
    gene_drop_rate = 0,
    work_dir = "./"
    )
```

##### （2）本地运行使用

```julia
#将RankCompV3包从github上clone到本地
git clone https://github.com/Yanjj1/RankCompV3.git
#在julia中加载RankCompV3包
#path表示git clone https://github.com/Yanjj1/RankCompV3.git到的路径
include("path/RankCompV3/src/RankCompV3.jl")
using Main.RankCompV3
#runing
reoa(fn_expr = "expr.txt",
        fn_metadata = "metadata.txt";
    expr_threshold = 3,
          pval_reo = 0.01,
     pval_sign_reo = 1.00,
     padj_sign_reo = 0.05,
           hk_file = "HK_genes_info.tsv",
           hk_name = "ENSEMBL",
      ref_gene_num = 3000,
    no_use_housekeeping = 0,
           species = "human",
    cell_drop_rate = 0,
    gene_drop_rate = 0,
    work_dir = "./"
    )
```

#### 2）软件包.exe

```
RankCompV3/bin/RankCompV3 




```

### 5、输出文件说明

#### 1）文件

##### （1）若使用管家基因，则生成以下7个文件。反之，则不生成。

###### 		a.在每组样本中随机抽取3个样本绘制管家基因在样本中的表达分布图.pdf文件。

###### 	b.基因表达谱.tsv文件，行表示基因，列表示样本。其中，第一列为基因名，第二列为是否属于管家基因。

##### （2）每次迭代产生6个文件

###### 	a.本次迭代得到的差异基因结果.tsv文件。其中，共5列，每列分别为gene name, delta, sd_delta, p, FDR.

###### 	b.The distribution maps of sd_delta, Pval and Padj of the whole genes.

###### 	c.本次迭代管家基因的基因表达谱.tsv文件，行表示基因，列表示样本。其中，第一列为基因名，第二列为是否属于管家基因。

##### （3）差异基因稳定后的最终输出文件

###### 	a.The final calculation results (gene name, delta, sd_delta, p, FDR) without threshold screening.

###### 	b.Differential gene expression profile.

###### 	c.Differential gene expression heat map file.

#### 2）log日志

### 6、应用示例

#### 1）code

#### 1）julia中RankCompV3包

[Yanjj1/RankCompV3: Differential expression gene recognition algorithm based on REOs. (github.com)](https://github.com/Yanjj1/RankCompV3)

##### （1）直接使用

```julia
#RankCompV3包下载方式详见4、工具使用方式
using RankCompV3
#包自带示例数据。使用默认参数，若需修改，可直接进行参数添加。
reoa_test()
#或本地文件使用。使用默认参数若需修改，可直接进行参数添加。
reoa("/public/yanj/data/fn_expr.txt",
	"/public/yanj/data/fn_metadata.txt"
)
```

##### （2）本地运行使用

```julia
#RankCompV3包下载方式详见4、工具使用方式
using Main.RankCompV3
#包自带示例数据。使用默认参数，若需修改，可直接进行参数添加。
reoa_test()
#或本地文件使用。使用默认参数若需修改，可直接进行参数添加。
reoa("/public/yanj/data/fn_expr.txt",
	"/public/yanj/data/fn_metadata.txt"
)
```

#### 2）软件包.exe

```
RankCompV3/bin/RankCompV3 




```



#### 2）输入文件

##### （1）表达谱文件fn_expr.txt

![image-20221008161103190](E:\资料\typora\图像\image-20221008161103190.png)

##### （2）metadata文件

![image-20221008161129940](E:\资料\typora\图像\image-20221008161129940.png)

##### （3）管家基因文件（内置，同时支持重新提供）

![image-20221007221716873](E:\资料\typora\图像\image-20221007221716873.png)

#### 3）输出文件

##### （1） 若使用管家基因，则生成以下7个文件。反之，则不生成。

###### 	a.在每组样本中随机抽取3个样本绘制管家基因在样本中的表达分布图.pdf文件。

![image-20221008160926863](E:\资料\typora\图像\image-20221008160926863.png)

###### 文件内容：

<center class="half">
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_hk_nonhk_gene_sample2.pdf" width="33%" height="250" >
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_hk_nonhk_gene_sample6.pdf" width="33%" height="250" >
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_hk_nonhk_gene_sample1.pdf" width="33%" height="250" >
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_hk_nonhk_gene_sample5.pdf" width="33%" height="250" >
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_hk_nonhk_gene_sample4.pdf" width="33%" height="250" >
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_hk_nonhk_gene_sample8.pdf" width="33%" height="250" >
</center>

###### b.基因表达谱.tsv文件，行表示基因，列表示样本。其中，第一列为基因名，第二列为是否属于管家基因。

![image-20221008161442156](E:\资料\typora\图像\image-20221008161442156.png)

###### 文件内容：

![image-20221008164912213](E:\资料\typora\图像\image-20221008164912213.png)

##### （2）每次迭代产生6个文件（以下仅对第0次迭代结果进行展示）

###### a.本次迭代得到的差异基因结果.tsv文件。其中，共5列，每列分别为gene name, delta, sd_delta, p, FDR.

![image-20221008172208671](E:\资料\typora\图像\image-20221008172208671.png)

###### 文件内容：

![image-20221008171910717](E:\资料\typora\图像\image-20221008171910717.png)

###### b.The distribution maps of sd_delta, Pval and Padj of the whole genes.
![image-20221008172449726](E:\资料\typora\图像\image-20221008172449726.png)

###### 文件内容：

<center class="half">
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_0_delta_graph.pdf" width="49%" height="430" >
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_0_sd_delta_graph.pdf" width="49%" height="430" >
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_0_Pval_graph.pdf" width="49%" height="430" >
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_0_Padj_graph.pdf" width="49%" height="430" >
</center>
###### c.本次迭代管家基因的基因表达谱.tsv文件，行表示基因，列表示样本。其中，第一列为基因名，第二列为是否属于管家基因。

![image-20221008172333440](E:\资料\typora\图像\image-20221008172333440.png)

###### 文件内容：

![image-20221008171949546](E:\资料\typora\图像\image-20221008171949546.png)

##### （3）差异基因稳定后的最终输出文件

###### a.The final calculation results (gene name, delta, sd_delta, p, FDR) without threshold screening.

![image-20221008193800156](E:\资料\typora\图像\image-20221008193800156.png)

###### 文件内容：

![image-20221008193740965](E:\资料\typora\图像\image-20221008193740965.png)

###### b.Differential gene expression profile.

![image-20221008193530085](E:\资料\typora\图像\image-20221008193530085.png)

###### 文件内容：

![image-20221008194050453](E:\资料\typora\图像\image-20221008194050453.png)

###### c.Differential gene expression heat map file.

![image-20221008194229822](E:\资料\typora\图像\image-20221008194229822.png)

###### 文件内容：
<center class="half">
    <embed id="pdfPlayer" src="../RankCompV3-test-data-output/RankCompV3_test_data_output/fn_expr_deg_exp_graph.pdf" width="100%" height="250" >
</center>















3）shell中运行.jl脚本

```shell
#注意，使用此方法，依赖包需要自行安装
julia test.jl --fn_expr "/home/yanj/jupyter_work/McCullagh/outcome_complete/data_outcome/rankcompV3_data_test/GSE41328_exp.txt" --fn_metadata "/home/yanj/jupyter_work/McCullagh/outcome_complete/data_outcome/rankcompV3_data_test/metadata.txt" --expr_threshold 3 --pval_reo 0.01 --pval_sign_reo 1.00 --padj_sign_reo 0.05 --hk_file "/home/yanj/jupyter_work/McCullagh/outcome_complete/HK_genes_info.tsv" --hk_name "ENSEMBL" --ref_gene_num 3000 --no_use_housekeeping 0 --species "human" --cell_drop_rate 0 --gene_drop_rate 0 1>test.log 2>test.err
```



```
INFO: The species information you selected for input data is: human
INFO: Successfully read in the expression matrix, including gene names (1st column),  with a size of: (43346, 9)
INFO: There were 1 samples with the number of detected genes (non-0 value) less than 0, and the samples were removed.
INFO: The number of cells detected in the presence of 7526 gene (value non-0) was less than 0 in the ctrl and/or treat groups, and gene was removed.
INFO: Size of the matrices after filtering non-expressed genes
INFO: Size of the control group: (14633, 4)
INFO: Size of the treatment group: (14633, 4)
INFO: Successfully read in the house-keeping genes, with a size of: (16, 1)
The number of housekeeping genes was 6, and the number of non-housekeeping genes was 14627. The expression profile (gene name, housekeeping gene or nsample) was saved into GSE82158_shielded_chimera_norm_counts_hk_nonhk_gene.tsv.
┌ Warning: RCall.jl: Warning: Removed 10378 rows containing non-finite values (stat_density).
└ @ RCall ~/.julia/packages/RCall/6kphM/src/io.jl:172
┌ Warning: RCall.jl: Warning: Removed 10408 rows containing non-finite values (stat_density).
└ @ RCall ~/.julia/packages/RCall/6kphM/src/io.jl:172
┌ Warning: RCall.jl: Warning: Removed 10344 rows containing non-finite values (stat_density).
└ @ RCall ~/.julia/packages/RCall/6kphM/src/io.jl:172
┌ Warning: RCall.jl: Warning: Removed 14588 rows containing non-finite values (stat_density).
└ @ RCall ~/.julia/packages/RCall/6kphM/src/io.jl:172
┌ Warning: RCall.jl: Warning: Removed 10341 rows containing non-finite values (stat_density).
└ @ RCall ~/.julia/packages/RCall/6kphM/src/io.jl:172
┌ Warning: RCall.jl: Warning: Removed 10232 rows containing non-finite values (stat_density).
└ @ RCall ~/.julia/packages/RCall/6kphM/src/io.jl:172
WARN: even if all samples have the identical REOs, it still cannot reach the required significance, 0.01.
      The max. significance is 0.125.
WARN: even if all samples have the identical REOs, it still cannot reach the required significance, 0.01.
      The max. significance is 0.125.
INFO: Minimum size for significantly stable REO in the control group: 4
INFO: Minimum size for significantly stable REO in the treatment group: 4
INFO: Number of threads, 1
 19.403717 seconds (45.71 M allocations: 2.215 GiB, 5.81% gc time, 87.02% compilation time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_0_result.tsv exists and will be overwritten.
INFO: N of 14311 genes in the 0th calculation is a singular matrix, and there are 0 differential genes.Save the results (gene name, delta, sd_delta,as GSE82158_shielded_chimera_norm_counts_iteration_0_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 0th iteration were drawn. Save as GSE82158_shielded_corm_counts_0_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_0_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_0_Padj_graph.pdf
INFO: the number of non-managed genes in the 0th iteration is 0, the expression profile (gene name, whether to manage gene, sample) has been saved as8_shielded_chimera_norm_counts_hk_nonhk_gene_0.tsv.
INFO: Number of threads, 1
^[[C^[[C^[[A300.944332 seconds (4.71 G allocations: 280.889 GiB, 20.84% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_1_result.tsv exists and will be overwritten.
INFO: N of 8 genes in the 1th calculation is a singular matrix, and there are 1808 differential genes.Save the results (gene name, delta, sd_delta, ps GSE82158_shielded_chimera_norm_counts_iteration_1_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 1th iteration were drawn. Save as GSE82158_shielded_corm_counts_1_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_1_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_1_Padj_graph.pdf
INFO: the number of non-managed genes in the 1th iteration is 1808, the expression profile (gene name, whether to manage gene, sample) has been saved2158_shielded_chimera_norm_counts_hk_nonhk_gene_1.tsv.
INFO: Number of threads, 1
267.002789 seconds (4.13 G allocations: 246.209 GiB, 21.27% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_2_result.tsv exists and will be overwritten.
INFO: N of 9 genes in the 2th calculation is a singular matrix, and there are 2560 differential genes.Save the results (gene name, delta, sd_delta, ps GSE82158_shielded_chimera_norm_counts_iteration_2_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 2th iteration were drawn. Save as GSE82158_shielded_corm_counts_2_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_2_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_2_Padj_graph.pdf
INFO: the number of non-managed genes in the 2th iteration is 2560, the expression profile (gene name, whether to manage gene, sample) has been saved2158_shielded_chimera_norm_counts_hk_nonhk_gene_2.tsv.
INFO: Number of threads, 1
249.021982 seconds (3.89 G allocations: 231.784 GiB, 21.30% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_3_result.tsv exists and will be overwritten.
INFO: N of 11 genes in the 3th calculation is a singular matrix, and there are 2805 differential genes.Save the results (gene name, delta, sd_delta,as GSE82158_shielded_chimera_norm_counts_iteration_3_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 3th iteration were drawn. Save as GSE82158_shielded_corm_counts_3_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_3_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_3_Padj_graph.pdf
INFO: the number of non-managed genes in the 3th iteration is 2805, the expression profile (gene name, whether to manage gene, sample) has been saved2158_shielded_chimera_norm_counts_hk_nonhk_gene_3.tsv.
INFO: Number of threads, 1
243.547466 seconds (3.81 G allocations: 227.087 GiB, 21.12% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_4_result.tsv exists and will be overwritten.
INFO: N of 10 genes in the 4th calculation is a singular matrix, and there are 2906 differential genes.Save the results (gene name, delta, sd_delta,as GSE82158_shielded_chimera_norm_counts_iteration_4_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 4th iteration were drawn. Save as GSE82158_shielded_corm_counts_4_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_4_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_4_Padj_graph.pdf
INFO: the number of non-managed genes in the 4th iteration is 2906, the expression profile (gene name, whether to manage gene, sample) has been saved2158_shielded_chimera_norm_counts_hk_nonhk_gene_4.tsv.
INFO: Number of threads, 1
241.184288 seconds (3.78 G allocations: 225.149 GiB, 20.95% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_5_result.tsv exists and will be overwritten.
INFO: N of 11 genes in the 5th calculation is a singular matrix, and there are 2949 differential genes.Save the results (gene name, delta, sd_delta,as GSE82158_shielded_chimera_norm_counts_iteration_5_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 5th iteration were drawn. Save as GSE82158_shielded_corm_counts_5_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_5_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_5_Padj_graph.pdf
INFO: the number of non-managed genes in the 5th iteration is 2949, the expression profile (gene name, whether to manage gene, sample) has been saved2158_shielded_chimera_norm_counts_hk_nonhk_gene_5.tsv.
INFO: Number of threads, 1
240.864387 seconds (3.76 G allocations: 224.325 GiB, 20.86% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_6_result.tsv exists and will be overwritten.
INFO: N of 9 genes in the 6th calculation is a singular matrix, and there are 2953 differential genes.Save the results (gene name, delta, sd_delta, ps GSE82158_shielded_chimera_norm_counts_iteration_6_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 6th iteration were drawn. Save as GSE82158_shielded_corm_counts_6_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_6_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_6_Padj_graph.pdf
INFO: the number of non-managed genes in the 6th iteration is 2953, the expression profile (gene name, whether to manage gene, sample) has been saved2158_shielded_chimera_norm_counts_hk_nonhk_gene_6.tsv.
INFO: Number of threads, 1
240.968264 seconds (3.76 G allocations: 224.247 GiB, 20.68% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_7_result.tsv exists and will be overwritten.
INFO: N of 10 genes in the 7th calculation is a singular matrix, and there are 2972 differential genes.Save the results (gene name, delta, sd_delta,as GSE82158_shielded_chimera_norm_counts_iteration_7_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 7th iteration were drawn. Save as GSE82158_shielded_corm_counts_7_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_7_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_7_Padj_graph.pdf
INFO: the number of non-managed genes in the 7th iteration is 2972, the expression profile (gene name, whether to manage gene, sample) has been saved2158_shielded_chimera_norm_counts_hk_nonhk_gene_7.tsv.
INFO: Number of threads, 1
241.266749 seconds (3.75 G allocations: 223.884 GiB, 20.65% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_8_result.tsv exists and will be overwritten.
INFO: N of 9 genes in the 8th calculation is a singular matrix, and there are 2983 differential genes.Save the results (gene name, delta, sd_delta, ps GSE82158_shielded_chimera_norm_counts_iteration_8_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 8th iteration were drawn. Save as GSE82158_shielded_corm_counts_8_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_8_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_8_Padj_graph.pdf
INFO: the number of non-managed genes in the 8th iteration is 2983, the expression profile (gene name, whether to manage gene, sample) has been saved2158_shielded_chimera_norm_counts_hk_nonhk_gene_8.tsv.
INFO: Number of threads, 1
241.339506 seconds (3.75 G allocations: 223.674 GiB, 20.54% gc time)
WARN: file GSE82158_shielded_chimera_norm_counts_iteration_9_result.tsv exists and will be overwritten.
INFO: N of 10 genes in the 9th calculation is a singular matrix, and there are 2983 differential genes.Save the results (gene name, delta, sd_delta,as GSE82158_shielded_chimera_norm_counts_iteration_9_result.tsv file except for the genes of singular matrix.
INFO: The distribution maps of sd_delta, Pval and Padj of the whole genes in the results of the 9th iteration were drawn. Save as GSE82158_shielded_corm_counts_9_sd_delta_graph.pdf、GSE82158_shielded_chimera_norm_counts_9_Pval_graph.pdf、和GSE82158_shielded_chimera_norm_counts_9_Padj_graph.pdf
INFO: After 9 iterations, the obtained list of differential genes (gene name, delta, sd_delta, p, FDR) is saved as GSE82158_shielded_chimera_norm_coult.tsv file, Save differential gene expression profile as GSE82158_shielded_chimera_norm_counts_deg_exp.tsv.
INFO: The final calculation results (gene name, delta, sd_delta, p, FDR) without threshold screening are GSE82158_shielded_chimera_norm_counts_iteratsult.tsv files.
INFO: The expression heat map of DEGs was drawn and saved as GSE82158_shielded_chimera_norm_counts_deg_exp_graph.pdf
2983×5 DataFrame
  Row │ names               delta     sd_delta  Pval          Padj
      │ String              Float64   Float64   Float64       Float64
──────┼────────────────────────────────────────────────────────────────────
    1 │ ENSMUSG00000000028  0.362627   6.31539  1.34738e-10   8.78051e-10
    2 │ ENSMUSG00000000093  0.443033  11.4647   9.92258e-31   9.53572e-30
    3 │ ENSMUSG00000000125  0.483575   9.0957   4.69916e-20   3.75342e-19
    4 │ ENSMUSG00000000127  0.214791   3.84485  6.03138e-5    0.000329449
    5 │ ENSMUSG00000000282  1.17752   17.3767   6.19136e-68   9.46643e-67
    6 │ ENSMUSG00000000301  0.646007  13.0987   1.67483e-39   1.86026e-38
    7 │ ENSMUSG00000000384  0.34975    6.36951  9.48195e-11   6.21793e-10
    8 │ ENSMUSG00000000402  0.108144   2.39683  0.0082687     0.0407222
    9 │ ENSMUSG00000000531  0.474683   8.68839  1.83812e-18   1.42611e-17
   10 │ ENSMUSG00000000579  1.17545   18.6163   1.18496e-77   2.02106e-76
   11 │ ENSMUSG00000000823  0.390372   6.09944  5.32192e-10   3.41486e-9
   12 │ ENSMUSG00000000902  0.594111   8.57752  4.84698e-18   3.72298e-17
   13 │ ENSMUSG00000000916  0.152703   2.71585  0.00330527    0.0167058
   14 │ ENSMUSG00000000934  0.371489   6.3261   1.25719e-10   8.2092e-10
   15 │ ENSMUSG00000000958  0.539068   8.19507  1.25225e-16   9.38676e-16
   16 │ ENSMUSG00000001029  2.55301   46.0771   0.0           0.0
  ⋮   │         ⋮              ⋮         ⋮           ⋮             ⋮
 2968 │ ENSMUSG00000103380  0.577827   9.8845   2.4299e-23    2.0674e-22
 2969 │ ENSMUSG00000103400  0.280895   5.19531  1.02192e-7    6.1496e-7
 2970 │ ENSMUSG00000103475  1.45721   31.1013   1.1564e-212   5.26576e-211
 2971 │ ENSMUSG00000103547  0.217257   4.77889  8.81339e-7    5.15095e-6
 2972 │ ENSMUSG00000103570  0.195213   3.88629  5.08934e-5    0.000278618
 2973 │ ENSMUSG00000103622  0.855458  17.9155   4.46373e-72   7.15421e-71
 2974 │ ENSMUSG00000103630  0.357111   8.77358  8.65405e-19   6.75367e-18
 2975 │ ENSMUSG00000103901  1.43502   30.3882   3.92719e-203  1.69834e-201
 2976 │ ENSMUSG00000103922  0.157237   3.34962  0.000404614   0.00213126
 2977 │ ENSMUSG00000104350  0.231987   5.33236  4.84733e-8    2.94119e-7
 2978 │ ENSMUSG00000104367  0.998009  21.1783   7.5698e-100   1.58521e-98
 2979 │ ENSMUSG00000104377  0.618824  11.5757   2.73503e-31   2.65634e-30
 2980 │ ENSMUSG00000104444  0.331103   6.42655  6.52645e-11   4.29911e-10
 2981 │ ENSMUSG00000104453  1.17714   19.2817   3.82553e-83   6.82756e-82
 2982 │ ENSMUSG00000104467  1.02923   22.8895   2.95818e-116  6.94056e-115
 2983 │ ENSMUSG00000104523  1.70206   24.8618   9.62817e-137  2.68578e-135
                                                          2951 rows omitted


```





人和小鼠物种，高通量测序数据（RNA-seq）、单细胞测序数据（scRNA-seq）、CEL-seq测序数据和蛋白质组数据分析。

edgeR是一个用于数字基因表达数据差异表达分析的包，即DNA测序技术产生的计数数据。它特别设计用于RNA-Seq或SAGE数据的差异表达分析，或ChIP-Seq数据的差异标记分析。



该差异基因识别工具基于julia语言开发，具体使用如下。



```julia
julia test.jl --fn_expr "/home/yanj/jupyter_work/McCullagh/outcome_complete/data_outcome/rankcompV3_data_test/GSE41328_exp.txt" --fn_metadata "/home/yanj/jupyter_work/McCullagh/outcome_complete/data_outcome/rankcompV3_data_test/metadata.txt" --hk_file "/home/yanj/jupyter_work/McCullagh/outcome_complete/HK_genes_info.tsv" 1>test.log 2>test.err
```

适用于高通量测序数据、scRNA-seq数据和蛋白质组数据等研究。

算法基于相对表达秩序关系REO设计，具体的算法原理可参看。。。



```julia
julia test.jl --fn_expr "/home/yanj/jupyter_work/McCullagh/outcome_complete/data_outcome/rankcompV3_data_test/GSE41328_exp.txt" --fn_metadata "/home/yanj/jupyter_work/McCullagh/outcome_complete/data_outcome/rankcompV3_data_test/metadata.txt" --expr_threshold 3 --pval_reo 0.01 --pval_sign_reo 1.00 --padj_sign_reo 0.05 --hk_file "/home/yanj/jupyter_work/McCullagh/outcome_complete/HK_genes_info.tsv" --hk_name "ENSEMBL" --ref_gene_num 3000 --no_use_housekeeping 0 --species "human" --cell_drop_rate 0 --gene_drop_rate 0 1>test.log 2>test.err
```







