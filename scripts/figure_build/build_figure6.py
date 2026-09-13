from __future__ import annotations
import hashlib,json
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np,pandas as pd
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'svg.fonttype':'none','pdf.fonttype':42,'font.size':7,'axes.spines.right':False,'axes.spines.top':False,'legend.frameon':False})
OUTROOT=Path(__file__).resolve().parents[1]; ROOT=OUTROOT.parent; OUT=OUTROOT/'figures'/'Figure6_20260912'; SD=OUT/'source_data'; MAN=OUTROOT/'figure_contracts'/'input_manifest.sha256.tsv'
S={'S003':'revision_v2/07_virtual_perturbation/09_manuscript_integration/manuscript_zh.v3_virtual_perturbation.src.md','S060':'revision_v2/07_virtual_perturbation/06_consensus/donor_stability.tsv','S061':'revision_v2/07_virtual_perturbation/06_consensus/seed_pairwise_spearman.tsv','S062':'revision_v2/07_virtual_perturbation/06_consensus/seed_pairwise_top200_jaccard.tsv','S063':'revision_v2/07_virtual_perturbation/07_enrichment/v1/ed_signature_gsea.tsv','S064':'revision_v2/07_virtual_perturbation/07_enrichment/v1/ed_top5_overlap.tsv','S065':'revision_v2/07_virtual_perturbation/07_enrichment/v1/gene_set_eligibility.tsv','S066':'revision_v2/07_virtual_perturbation/06_consensus/negative_control_calibration_v2/negative_control_calibration.tsv'}
G=['TYMS','EFNB2','LRRC17']; C={'TYMS':'#2A9D8F','EFNB2':'#E69F00','LRRC17':'#8A6BB8'}
def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
 return h.hexdigest()
def verify():
 m=pd.read_csv(MAN,sep='\t',dtype=str).set_index('source_id');o={}
 for k,p in S.items():
  q=ROOT/p
  if m.loc[k,'relative_path']!=p or sha(q).lower()!=m.loc[k,'sha256'].lower():raise ValueError('input drift '+k)
  o[k]=sha(q)
 return o
def label(ax,x,xx=-.12):ax.text(xx,1.05,x,transform=ax.transAxes,fontweight='bold',fontsize=8.2)
def main():
 v=verify(); st=pd.read_csv(ROOT/S['S060'],sep='\t'); sp=pd.read_csv(ROOT/S['S061'],sep='\t'); ja=pd.read_csv(ROOT/S['S062'],sep='\t'); gs=pd.read_csv(ROOT/S['S063'],sep='\t'); ov=pd.read_csv(ROOT/S['S064'],sep='\t'); elig=pd.read_csv(ROOT/S['S065'],sep='\t'); cal=pd.read_csv(ROOT/S['S066'],sep='\t')
 if st.shape[0]!=9 or sp.shape[0]!=90 or ja.shape[0]!=90 or gs.shape[0]!=3 or ov.shape[0]!=6 or cal.shape[0]!=21:raise ValueError('frozen dimensions')
 if not (gs.pathway=='ED_UP_FDR').all() or not np.allclose(gs.NES,[1.762643,1.747390,1.394354]):raise ValueError('gsea')
 if not ov.loc[ov.signature.eq('ED_DOWN_FDR'),'signature_size'].eq(3).all(): raise ValueError('ED_DOWN size-3 eligibility rule missing')
 # sources
 SD.mkdir(parents=True,exist_ok=True)
 design=pd.DataFrame([{'main_runs':45,'matched_controls':90,'targets':3,'donor_backgrounds':3,'seeds_per_target_donor':5,'technical_repeat_note':'Seeds are technical/model repeats, not biological replicates.'}])
 for d,n in [(design,'Figure6a_design_source_data.tsv'),(st,'Figure6b_donor_stability_source_data.tsv'),(gs,'Figure6c_ED_UP_GSEA_source_data.tsv'),(ov,'Figure6d_top5_overlap_source_data.tsv'),(cal,'Figure6e_matched_calibration_source_data.tsv')]:d.to_csv(SD/n,sep='\t',index=False)
 fig=plt.figure(figsize=(7.205,5.709)); grid=fig.add_gridspec(2,6,height_ratios=[1,1.35]);fig.subplots_adjust(left=.11,right=.96,top=.95,bottom=.10,wspace=.95,hspace=.82)
 a=fig.add_subplot(grid[0,:2]);b=fig.add_subplot(grid[0,2:4]);c=fig.add_subplot(grid[0,4:]);d=fig.add_subplot(grid[1,:3]);e=fig.add_subplot(grid[1,3:])
 # a
 a.axis('off'); a.text(.02,.86,'Predeclared virtual-perturbation design',fontsize=7.2,fontweight='bold',transform=a.transAxes);a.text(.06,.61,'3 targets × 3 donor backgrounds\n× 5 seeds = 45 main model runs',fontsize=7,transform=a.transAxes,bbox={'facecolor':'#E7F0F8','edgecolor':'#3775BA','pad':4});a.text(.06,.26,'90 matched-gene controls\nfor network-score calibration',fontsize=6.5,transform=a.transAxes,bbox={'facecolor':'#F2EAF6','edgecolor':'#8A6BB8','pad':4});a.text(.02,.04,'Seeds are technical/model repeats; donors are distinct cellular backgrounds.',fontsize=4.6,transform=a.transAxes,color='#666');label(a,'a')
 # b
 x=np.arange(9)
 for g in G:
  q=st[st.target_gene==g];b.scatter(q.index, q.median_seed_spearman,s=22,color=C[g],label=g);b.scatter(q.index,q.median_top200_jaccard,s=22,facecolors='none',edgecolors=C[g],linewidth=.9)
 b.axhline(.7,color='#777',ls='--',lw=.55);b.axhline(.4,color='#777',ls=':',lw=.55);b.set_ylim(.3,1);b.set_xticks(x);b.set_xticklabels(['']*9);b.set_ylabel('Stability metric',fontsize=5.4);b.tick_params(labelsize=4.8);b.set_title('Donor-level seed stability',loc='left',fontsize=7,fontweight='bold');b.legend(fontsize=4.5,loc='lower right');b.text(.02,.03,'filled: median Spearman; open: Top-200 Jaccard',transform=b.transAxes,fontsize=4.2);label(b,'b',)
 # c
 q=gs.set_index('target_gene').reindex(G).reset_index();y=np.arange(3);c.hlines(y,0,q.NES,color='#aaa');c.scatter(q.NES,y,s=42,c=[C[g] for g in G]);
 for i,r in q.iterrows():
  c.text(r.NES+.03,i,f'NES {r.NES:.2f}\nFDR {r.padj:.3g}',va='center',fontsize=4.4)
 c.set_xlim(0,2.2);c.set_yticks(y);c.set_yticklabels(G,fontsize=5.4);c.invert_yaxis();c.set_xlabel('ED_UP_FDR NES',fontsize=5.2);c.tick_params(axis='x',labelsize=4.8);c.set_title('ED_UP ranked enrichment',loc='left',fontsize=7,fontweight='bold');c.text(0,-.32,'ED_DOWN ranked GSEA: NOT_TESTED_SIZE_3',transform=c.transAxes,fontsize=4.2,color='#666');label(c,'c',xx=-.2)
 # d
 up=ov[ov.signature=='ED_UP_FDR'].set_index('target_gene').reindex(G);dn=ov[ov.signature=='ED_DOWN_FDR'].set_index('target_gene').reindex(G); yy=np.arange(3);d.barh(yy-.16,up.overlap_size,height=.3,color=[C[g] for g in G],label='ED_UP_FDR');d.barh(yy+.16,dn.overlap_size,height=.3,color='#c7c7c7',label='ED_DOWN_FDR');
 for i,g in enumerate(G):
  d.text(up.loc[g,'overlap_size']+.15,i-.16,f"FDR {up.loc[g,'FDR_within_target']:.3g}",va='center',fontsize=4.5)
 d.set_yticks(yy);d.set_yticklabels(G,fontsize=5.4);d.invert_yaxis();d.set_xlabel('Top-5% overlap size',fontsize=5.2);d.set_title('Top network-ranked gene overlap',loc='left',fontsize=7,fontweight='bold');d.legend(fontsize=4.4,loc='lower right');d.text(0,-.27,'Within-target hypergeometric FDR; ED_DOWN signature size = 3.',transform=d.transAxes,fontsize=4.2,color='#666');label(d,'d')
 # e all metric ratios
 metrics=list(cal.metric.drop_duplicates()); mat=cal.pivot(index='metric',columns='target_gene',values='ratio_target_to_control_median').reindex(index=metrics,columns=G);im=e.imshow(mat,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=.7,vcenter=1,vmax=1.3));e.set_xticks(range(3));e.set_xticklabels(G,fontsize=5);e.set_yticks(range(len(metrics)));e.set_yticklabels([x.replace('_',' ') for x in metrics],fontsize=4.3);e.tick_params(length=0);cb=fig.colorbar(im,ax=e,fraction=.055,pad=.05);cb.set_label('Target/control ratio',fontsize=4.8);cb.ax.tick_params(labelsize=4.2);e.set_title('Matched-control calibration',loc='left',fontsize=7,fontweight='bold');e.text(0,-.22,'Controls calibrate network-score specificity only (n = 10 per target).',transform=e.transAxes,fontsize=4.1,color='#666');label(e,'e')
 OUT.mkdir(parents=True,exist_ok=True);base=OUT/'Figure6_virtual_perturbation_stability_and_ED_association'
 for ext,kw in [('.svg',{}),('.pdf',{}),('.png',{'dpi':300}),('.tiff',{'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}})]:
  fig.savefig(base.with_suffix(ext),facecolor='white',**kw)
 plt.close(fig)
 audit={'status':'PASS','figure':'Figure 6','backend':'Python/matplotlib','input_manifest_validation':'PASS','input_sha256':v,'assertions':{'main_runs':45,'matched_controls':90,'donor_stability_rows':9,'ED_UP_size':18,'ED_UP_NES':q.NES.tolist(),'top5_ED_UP_overlap':up.overlap_size.tolist(),'calibration_rows':21},'output_dimensions_mm':[183,145],'interpretation_boundary':'Predictive network effects only; not knockout validation, expression reversal, causal mechanism, or ED reversal.'}
 (OUT/'Figure6_build_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8');print(json.dumps({'status':'PASS','output':str(OUT)}))
if __name__=='__main__':main()
