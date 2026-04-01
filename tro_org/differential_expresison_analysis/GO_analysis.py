import gseapy as gp
import numpy as np
from matplotlib import pyplot as plt
import scanpy as sc
import seaborn as sns
import textwrap

data_res = {'Shibata': ['CGA', 'MALAT1', 'CSH1', 'RPS2', 'RPLP1', 'RPS18', 'RPS8', 'RPS19',
       'RPL41', 'RPS14', 'RPS12', 'RPS27', 'RPS28', 'RPL37', 'RPL10', 'RPL39',
       'RPS6', 'RPL34', 'RPS3', 'RPL13', 'RPL28', 'RPS4X', 'RPS24', 'RPL37A',
       'RPL32', 'RPL11', 'RPL8', 'RPL7A', 'RPS3A', 'RPL15', 'RPS29', 'RPL18',
       'RPL18A', 'RPS27A', 'RPLP2', 'RPL12', 'RPLP0', 'RPL19', 'RPL6', 'RPL36',
       'RPL3', 'RPS16', 'RPS11', 'RPL13A', 'RPS15', 'RPS15A', 'RPS13', 'RPL29',
       'RPS23', 'RPL9', 'RPL26', 'RPS7', 'RPL30', 'RPSA', 'RPL23A', 'RPL23',
       'RPL27A', 'RPL21', 'RPL14', 'RPL10A', 'RPS9', 'RPS5', 'RPL7', 'RPS25',
       'RPL5', 'RPL35A', 'RPL24', 'RPL35', 'RPS21', 'RPS20', 'RPL22', 'RPL27',
       'RPL4', 'RPS26', 'RPL38', 'RPL31', 'PAGE4', 'RPS27L', 'XAGE3', 'RPS10',
       'RPL36AL', 'AES', 'HSD3B1', 'RPL22L1', 'RPL36A', 'SEPT11', 'PLA2G16',
       'SEPT7', 'RPS4Y1', 'SEPT2', 'RPS19BP1', 'FAM96B', 'C19orf70', 'C8orf59',
       'MINOS1', 'CSH2', 'INSL4', 'RPL17', 'LGALS16', 'DIRC2'],
'Arutyunyan_PTO': ['CD70', 'SLC25A21', 'TMSB15B', 'NKD2', 'CCND2', 'LRRN1', 'SALL4',
       'LDOC1', 'ATP2A3', 'SLC6A12', 'CREB3L1', 'PACSIN1', 'BNC1', 'PLD5',
       'DHDH', 'GAL3ST1', 'WNT3A', 'PRPS2', 'EHD2', 'SOX9', 'CALB2', 'F2R',
       'HECW2', 'TFF1', 'HS3ST1', 'COL1A2', 'LINC00632', 'SLC9A3R2', 'CTH',
       'FAM189A2', 'SH3GL3', 'SNHG18', 'NKAIN4', 'CRIP1', 'CD74', 'CDA',
       'IQCD', 'LMOD1', 'E2F8', 'EPCAM', 'NRTN', 'DBN1', 'PRSS23', 'HIST1H1A',
       'CCBE1', 'SPEG', 'RTKN2', 'CGN', 'C6orf132', 'CLDN11', 'TEKT2', 'NOP2',
       'SMPDL3B', 'GATA6', 'BDKRB1', 'OSBPL6', 'GAS6', 'C11orf96', 'ABHD6',
       'CARD10', 'SCUBE3', 'TMEM45A', 'RASL10A', 'SPATA9', 'JUNB', 'RHOBTB3',
       'STEAP4', 'PRSS22', 'KISS1R', 'NAV2-AS2', 'MMP11', 'EGFR-AS1'],
'Arutyunyan_TSC': ['IGFBP3', 'SLC16A12', 'TLL1', 'AKR1B15', 'COL14A1', 'TMEM218', 'PDE10A',
       'DHRS4', 'DHRS4L2', 'CTSS', 'ACADL', 'RAB38', 'FBXO15', 'RAB27B',
       'SLC2A9', 'EPHA1-AS1', 'DEFB1', 'IFI16', 'C1R', 'TMEM150C', 'ITGB6',
       'SERP2', 'PRRT1', 'APLN', 'ABCA12', 'GCNT4', 'ANKRD37', 'ACSS1',
       'SGMS2', 'BARX2', 'PIK3AP1', 'MAB21L3', 'ARHGAP24', 'CRNDE', 'VILL',
       'APOL1', 'ANGPTL4', 'GHR', 'MAN1A1', 'OPRK1', 'BCL2', 'C9orf116',
       'ZNF175', 'CCDC68', 'CASP1', 'RBM20', 'MME', 'VAV3', 'RARRES1',
       'FBXO32', 'TMEM91', 'NRG4', 'CD40', 'SPTLC3', 'STEAP3', 'NRP2',
       'B3GNT3', 'TMEM140', 'ZNF808', 'ZNF525', 'L1CAM', 'LLGL1', 'DBNDD1',
       'APOBEC3C', 'TUBA1A', 'EEF1A2', 'RBP7', 'UCHL1', 'LINC01405',
       'DCAF4L2'],
'Shannon': ['RACK1', 'ATP5F1E', 'MYDGF', 'ATP5MG', 'SELENOW', 'ATP5MD', 'ATP5MC2',
       'BEX3', 'SEM1', 'ATP5ME', 'JPT1', 'TCEAL9', 'ATP5F1B', 'SNU13',
       'ATP5F1C', 'ATP5F1A', 'TRIR', 'ERO1A', 'AFDN', 'MESD', 'SF3B6',
       'ATP5PO', 'SELENOS', 'ATP5IF1', 'LHFPL6', 'SARAF', 'NSD3', 'TENT5A',
       'RTRAF', 'UCA1', 'ELOC', 'CD24', 'SELENOF', 'MINDY2', 'GRAMD2B',
       'STMP1', 'SINHCAF', 'ATP5PB', 'NBDY', 'ADGRG6', 'VSIR', 'OGA', 'NDUFB8',
       'HACD3', 'SELENOT', 'NOP53', 'JPT2', 'TUT4', 'NRDC', 'TMBIM4',
       'SELENOK', 'RSRP1', 'MRPL57', 'NUP58', 'NECTIN3', 'PPP4R3A', 'LSM8',
       'REX1BD', 'CRYBG1', 'CYTOR', 'SHTN1', 'MMP24OS', 'ERG28', 'ZNF518A',
       'WASHC4', 'SELENOH', 'TUT7', 'METTL26', 'PUM3', 'CCDC186', 'NDUFAF8',
       'RTF2', 'NECTIN2', 'PCNX4', 'TMEM263', 'KMT5B', 'MTREX', 'WAPL', 'UFD1',
       'SEC22B', 'BUD23', 'PPP4R3B', 'SMIM26', 'ECPAS', 'RAB5IF', 'CEMIP2',
       'PAXX', 'NECTIN4', 'AC022784.1', 'RTL8C', 'FAAP20', 'SPART', 'CEBPZOS',
       'MFSD14C', 'TCAF1', 'MIGA1', 'RIPOR2', 'PUDP', 'GCN1', 'THUMPD3-AS1']}


def gsea_analysis(data, ds):
       enr = gp.enrichr(gene_list=list(data[ds]),
                        gene_sets=['GO_Biological_Process_2025', 'KEGG_2021_Human'],
                        organism='human',
                        outdir=None)
       return enr

def plot_go_analysis(enr, key, go_plot_dir, top_n=15):
       df = enr.results.sort_values("Adjusted P-value").head(top_n)
       df["neglog10_p"] = -np.log10(df["Adjusted P-value"])
       df["Term_wrapped"] = df["Term"].apply(
           lambda x: "\n".join(textwrap.wrap(str(x), width=35))
       )

       plt.figure(figsize=(10, 10))
       sns.barplot(
           data=df,
           x="neglog10_p",
           y="Term_wrapped",
           color="#383a6b"
       )

       plt.xlabel("-log10(adjusted p-value)")
       plt.ylabel("")
       plt.title(f"GO enrichment: {key}")
       plt.subplots_adjust(left=0.45)
       plt.tight_layout()
       plt.savefig(f"{go_plot_dir}/{key}_go_analysis.png")
       plt.show()

def scanpy_markergenes(adata):
       sc.tl.rank_genes_groups(adata, groupby='sample', method='wilcoxon')
       sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, values_to_plot='logfoldchanges')

       samples = adata.obs["sample"].unique()

       marker_dict = {}

       for s in samples:
              df = sc.get.rank_genes_groups_df(adata, group=s)
              sig = (df[(df["pvals_adj"] < 0.01) & (df["logfoldchanges"] > 2)]
                     .sort_values("scores", ascending=False)
                     .head(200))
              marker_dict[s] = sig["names"].tolist()

       print(marker_dict)
       enr = gsea_analysis(marker_dict, "Shannon")
       plot_go_analysis(enr, "Shannon")


def go_analysis(marker_dict, go_plot_dir):
       # adata_file = "database/integrated_data/merged_integration_final_integrated.h5ad"
       # adata = sc.read_h5ad(adata_file)
       # scanpy_markergenes(adata)
       for key in marker_dict.keys():
              result = gsea_analysis(marker_dict, key)
              plot_go_analysis(result, key, go_plot_dir)

if __name__ == '__main__':
    go_analysis(data_res, "tro_org/differential_expresison_analysis/GO")
