
import pyranges as pr
import pandas as pd
import numpy as np


rare_alleles_cases = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/joep-bpd-vars-export2_hg38_cases_hpo.bed.gz")
rare_alleles_cases.columns = ['Chromosome', 'Start', 'End', 'Ref', 'Alt', 'AC']
rare_alleles_cases.AC = rare_alleles_cases.AC.astype(int)
rare_alleles_cases = rare_alleles_cases[['AC']]

rare_alleles_controls = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/WGS10K_rare_alleles/joep-bpd-vars-export2_hg38_controls_hpo.bed.gz")
rare_alleles_controls.columns = ['Chromosome', 'Start', 'End', 'Ref', 'Alt', 'AC']
rare_alleles_controls.AC = rare_alleles_controls.AC.astype(int)
rare_alleles_controls = rare_alleles_controls[['AC']]


MK_humpies = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/humpy_regions_hg38/MK_humpies_annotated_with_SCARS_clf_hg38.bed")
MK_humpies.columns = ['Chromosome', 'Start', 'End', 'median_SCARS_clf']

EB_humpies = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/humpy_regions_hg38/EB_humpies_annotated_with_SCARS_clf_hg38.bed")
EB_humpies.columns = ['Chromosome', 'Start', 'End', 'median_SCARS_clf']

MON_humpies = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/humpy_regions_hg38/MON_humpies_annotated_with_SCARS_clf_hg38.bed")
MON_humpies.columns = ['Chromosome', 'Start', 'End', 'median_SCARS_clf']

B_humpies = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/humpy_regions_hg38/B_humpies_annotated_with_SCARS_clf_hg38.bed")
B_humpies.columns = ['Chromosome', 'Start', 'End', 'median_SCARS_clf']

rCD4_humpies = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/humpy_regions_hg38/rCD4_humpies_annotated_with_SCARS_clf_hg38.bed")
rCD4_humpies.columns = ['Chromosome', 'Start', 'End', 'median_SCARS_clf']

aCD4_humpies = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/humpy_regions_hg38/aCD4_humpies_annotated_with_SCARS_clf_hg38.bed")
aCD4_humpies.columns = ['Chromosome', 'Start', 'End', 'median_SCARS_clf']


rare_alleles_cases.ix = np.arange(len(rare_alleles_cases))

df_case = pd.DataFrame(0, index=np.arange(len(rare_alleles_cases)), columns=['AC', 'MK_id', 'EB_id', 'MON_id', 'B_id', 'rCD4_id', 'aCD4_id', 'median_SCARS_clf', 'case_control'])

df_case.loc[rare_alleles_cases.join(MK_humpies).ix, 'MK_id'] = 1
df_case.loc[rare_alleles_cases.join(MK_humpies).ix, 'median_SCARS_clf'] = rare_alleles_cases.join(MK_humpies).median_SCARS_clf.values

df_case.loc[rare_alleles_cases.join(EB_humpies).ix, 'EB_id'] = 1
df_case.loc[rare_alleles_cases.join(EB_humpies).ix, 'median_SCARS_clf'] = rare_alleles_cases.join(EB_humpies).median_SCARS_clf.values

df_case.loc[rare_alleles_cases.join(MON_humpies).ix, 'MON_id'] = 1
df_case.loc[rare_alleles_cases.join(MON_humpies).ix, 'median_SCARS_clf'] = rare_alleles_cases.join(MON_humpies).median_SCARS_clf.values

df_case.loc[rare_alleles_cases.join(B_humpies).ix, 'B_id'] = 1
df_case.loc[rare_alleles_cases.join(B_humpies).ix, 'median_SCARS_clf'] = rare_alleles_cases.join(B_humpies).median_SCARS_clf.values

df_case.loc[rare_alleles_cases.join(rCD4_humpies).ix, 'rCD4_id'] = 1
df_case.loc[rare_alleles_cases.join(rCD4_humpies).ix, 'median_SCARS_clf'] = rare_alleles_cases.join(rCD4_humpies).median_SCARS_clf.values

df_case.loc[rare_alleles_cases.join(aCD4_humpies).ix, 'aCD4_id'] = 1
df_case.loc[rare_alleles_cases.join(aCD4_humpies).ix, 'median_SCARS_clf'] = rare_alleles_cases.join(aCD4_humpies).median_SCARS_clf.values

df_case = df_case[df_case.iloc[:,1:7].sum(axis=1) > 0]
df_case.AC = rare_alleles_cases.as_df().loc[df_case.index, 'AC']

df_case.case_control = 1




rare_alleles_controls.ix = np.arange(len(rare_alleles_controls))

df_control = pd.DataFrame(0, index=np.arange(len(rare_alleles_controls)), columns=['AC', 'MK_id', 'EB_id', 'MON_id', 'B_id', 'rCD4_id', 'aCD4_id', 'median_SCARS_clf', 'case_control'])

df_control.loc[rare_alleles_controls.join(MK_humpies).ix, 'MK_id'] = 1
df_control.loc[rare_alleles_controls.join(MK_humpies).ix, 'median_SCARS_clf'] = rare_alleles_controls.join(MK_humpies).median_SCARS_clf.values

df_control.loc[rare_alleles_controls.join(EB_humpies).ix, 'EB_id'] = 1
df_control.loc[rare_alleles_controls.join(EB_humpies).ix, 'median_SCARS_clf'] = rare_alleles_controls.join(EB_humpies).median_SCARS_clf.values

df_control.loc[rare_alleles_controls.join(MON_humpies).ix, 'MON_id'] = 1
df_control.loc[rare_alleles_controls.join(MON_humpies).ix, 'median_SCARS_clf'] = rare_alleles_controls.join(MON_humpies).median_SCARS_clf.values

df_control.loc[rare_alleles_controls.join(B_humpies).ix, 'B_id'] = 1
df_control.loc[rare_alleles_controls.join(B_humpies).ix, 'median_SCARS_clf'] = rare_alleles_controls.join(B_humpies).median_SCARS_clf.values

df_control.loc[rare_alleles_controls.join(rCD4_humpies).ix, 'rCD4_id'] = 1
df_control.loc[rare_alleles_controls.join(rCD4_humpies).ix, 'median_SCARS_clf'] = rare_alleles_controls.join(rCD4_humpies).median_SCARS_clf.values

df_control.loc[rare_alleles_controls.join(aCD4_humpies).ix, 'aCD4_id'] = 1
df_control.loc[rare_alleles_controls.join(aCD4_humpies).ix, 'median_SCARS_clf'] = rare_alleles_controls.join(aCD4_humpies).median_SCARS_clf.values

df_control = df_control[df_control.iloc[:,1:7].sum(axis=1) > 0]
df_control.AC = rare_alleles_controls.as_df().loc[df_control.index, 'AC']

df_control.case_control = 0



df = pd.concat([df_case, df_control]).reset_index(drop=True)
df = df.loc[df.index.repeat(df.AC)].drop(['AC'],axis=1).reset_index(drop=True)
df.to_csv("/home/jwt44/case_control_regression_data.tsv", sep='\t', index=False)


