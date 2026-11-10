"""Build a shareable PDF from the editable, page-structured research report."""
from pathlib import Path
import re
import json
import hashlib
import shutil
import zipfile
from xml.sax.saxutils import escape
from PIL import Image as PILImage
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,Image,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT,TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=ROOT/'output/pdf'; OUT.mkdir(parents=True,exist_ok=True)
PDF=OUT/'embryo_twin_research_report.pdf'
for name,file in [('Arial','arial.ttf'),('Arial-Bold','arialbd.ttf'),('Arial-Italic','ariali.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(Path('C:/Windows/Fonts')/file)))
pdfmetrics.registerFontFamily('Arial',normal='Arial',bold='Arial-Bold',italic='Arial-Italic',boldItalic='Arial-Bold')
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='BodyReport',fontName='Arial',fontSize=9.7,leading=13.0,spaceAfter=6.5,textColor=colors.black))
styles.add(ParagraphStyle(name='ReportTitle',fontName='Arial-Bold',fontSize=24,leading=29,spaceAfter=14,textColor=colors.black))
styles.add(ParagraphStyle(name='SectionTitle',fontName='Arial-Bold',fontSize=17,leading=21,spaceAfter=12,textColor=colors.black))
styles.add(ParagraphStyle(name='Subheading',fontName='Arial-Bold',fontSize=11.8,leading=15,spaceBefore=7,spaceAfter=6,textColor=colors.black,keepWithNext=True))
styles.add(ParagraphStyle(name='CaptionReport',fontName='Arial',fontSize=8.3,leading=11,spaceBefore=3,spaceAfter=7,textColor=colors.HexColor('#333333')))
styles.add(ParagraphStyle(name='CellReport',fontName='Arial',fontSize=8.5,leading=11.2,spaceAfter=0,textColor=colors.black))
styles.add(ParagraphStyle(name='CellHeader',parent=styles['CellReport'],fontName='Arial-Bold',textColor=colors.white))
styles.add(ParagraphStyle(name='ListReport',parent=styles['BodyReport'],leftIndent=13,firstLineIndent=-10,spaceAfter=5))


def inline(s):
    s=escape(s)
    return re.sub(r'(https?://[^\s]+)',lambda m:'<link href="'+m.group(1)+'" color="#174A75">'+m.group(1)+'</link>',s)


class Report(SimpleDocTemplate):
    def afterFlowable(self,flowable):
        if getattr(flowable,'report_heading',False):
            key='section-'+str(len(self.heading_pages))
            self.canv.bookmarkPage(key); self.canv.addOutlineEntry(flowable.getPlainText(),key,0,False)
            self.heading_pages.append(dict(title=flowable.getPlainText(),page=self.page))


def footer(canvas,doc):
    canvas.saveState(); canvas.setFont('Arial',8); canvas.setFillColor(colors.HexColor('#555555'))
    canvas.drawString(46,25,'Embryo Twin research progress | 8 November 2026')
    canvas.drawRightString(A4[0]-46,25,str(doc.page)); canvas.restoreState()


def make_table(lines,width):
    parsed=[[part.strip() for part in line.strip().strip('|').split('|')] for line in lines]
    data=[parsed[0]]+parsed[2:]
    n=len(data[0]); first=data[0][0]
    if n==2:
        widths=[.14,.86] if first=='Record' else [.28,.72]
    elif n==3 and first=='Record': widths=[.12,.44,.44]
    elif n==3: widths=[.2,.36,.44]
    elif n==4: widths=[.37,.21,.21,.21] if first=='Condition' else [.26,.24,.25,.25]
    elif n==5: widths=[.13,.2,.22,.2,.25]
    else: widths=[1/n]*n
    cells=[[Paragraph(inline(c),styles['CellHeader' if i==0 else 'CellReport']) for c in row] for i,row in enumerate(data)]
    table=Table(cells,colWidths=[width*x for x in widths],repeatRows=1,hAlign='LEFT')
    commands=[('BACKGROUND',(0,0),(-1,0),colors.HexColor('#243F55')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),
              ('GRID',(0,0),(-1,-1),.45,colors.HexColor('#D9D9D9')),('LEFTPADDING',(0,0),(-1,-1),6),
              ('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]
    for i in range(1,len(data)):
        if i%2==0: commands.append(('BACKGROUND',(0,i),(-1,i),colors.HexColor('#F3F6F8')))
    table.setStyle(TableStyle(commands)); return table


def build():
    doc=Report(str(PDF),pagesize=A4,rightMargin=46,leftMargin=46,topMargin=42,bottomMargin=43,
        title='Self organized branching in a cell based simulation',author='Embryo Twin project',subject='Research progress and supervisor discussion')
    doc.heading_pages=[]; story=[]; text=(HERE/'report.md').read_text(encoding='utf-8')
    pages=text.split('<!-- PAGE -->')
    for number,page in enumerate(pages):
        if number: story.append(PageBreak())
        lines=page.strip().splitlines(); i=0
        while i<len(lines):
            line=lines[i].strip()
            if not line: i+=1; continue
            if line.startswith('# '):
                p=Paragraph(inline(line[2:]),styles['ReportTitle' if number==0 else 'SectionTitle']); p.report_heading=True
                story.append(p); i+=1; continue
            if line.startswith('## '): story.append(Paragraph(inline(line[3:]),styles['Subheading'])); i+=1; continue
            if line.startswith('|'):
                table=[]
                while i<len(lines) and lines[i].strip().startswith('|'): table.append(lines[i]); i+=1
                story.extend([make_table(table,doc.width),Spacer(1,9)]); continue
            match=re.match(r'!\[(.*?)\]\((.*?)\)',line)
            if match:
                path=HERE/match.group(2); w,h=PILImage.open(path).size
                scale=min(doc.width/w,205/h)
                image=Image(str(path),width=w*scale,height=h*scale)
                story.extend([Spacer(1,5),image,Spacer(1,5)]); i+=1; continue
            if line.startswith('- ') or re.match(r'^\d+\. ',line):
                story.append(Paragraph(inline((' - '+line[2:]) if line.startswith('- ') else line),styles['ListReport'])); i+=1; continue
            para=[line]; i+=1
            while i<len(lines) and lines[i].strip() and not lines[i].strip().startswith(('#','|','![')):
                para.append(lines[i].strip()); i+=1
            paragraph=' '.join(para)
            story.append(Paragraph(inline(paragraph),styles['CaptionReport' if paragraph.startswith('Figure ') else 'BodyReport']))
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    reader=PdfReader(PDF)
    (HERE/'layout.json').write_text(json.dumps(dict(pages=len(reader.pages),planned_sections=len(pages),heading_pages=doc.heading_pages),indent=2)+'\n')
    assert len(reader.pages)==len(pages),f'Overflow: {len(reader.pages)} pages for {len(pages)} planned pages'
    assert [h['page'] for h in doc.heading_pages]==list(range(1,len(pages)+1)),doc.heading_pages
    print('Created',PDF,'pages',len(reader.pages))


def evidence():
    # Bundle human-readable notebooks and selected summaries; retain raw NPZs in the repository.
    exp=ROOT/'research/experiments'
    notebooks=sorted(exp.glob('*/README.md'))
    files=[ROOT/'tools/BRANCHING_BASELINE.md',ROOT/'research/README.md',ROOT/'source/apps/embryo_twin.digital.kit']
    files+=sorted((ROOT/'source/extensions/embryo_twin.sim/embryo_twin/sim').glob('*.py'))
    files+=sorted((ROOT/'tools').glob('branching_*.py'))
    files+=sorted((ROOT/'tools').glob('test_branching_*.py'))
    for readme in notebooks:
        files.append(readme)
        folder=readme.parent
        files+=sorted(folder.glob('*.py'))
        for name in ['final_summary.csv','condition_summary.csv','resolution_differences.csv','validation.json','settings.json','sha256.json']:
            f=folder/name
            if f.exists(): files.append(f)
        if (folder/'analysis').exists():
            for name in ['final_summary.csv','window_summary.csv','final_summary.json']:
                f=folder/'analysis'/name
                if f.exists(): files.append(f)
    # Include the small raw inputs needed to reconstruct the four report figures.
    files.append(exp/'2026-11-07-pair-specific-adhesion/timecourses.csv')
    for topology in ['initial_sheet','contracted_E015']:
        files.append(exp/f'2026-11-07-static-tip-selection/seed42_{topology}_E015_dt0.025.npz')
    for timing in ['immediate','prepatterned']:
        files.append(exp/f'2026-11-08-junction-sprouting/seed42_{timing}_steps32_both.npz')
    files=list(dict.fromkeys(files))
    manifest={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    manifest.update({str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [HERE/'report.md',PDF,*sorted((HERE/'figures').glob('*.png'))]})
    (HERE/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    bundle=OUT/'embryo_twin_research_evidence.zip'
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
        for f in files: z.write(f,str(f.relative_to(ROOT)))
        z.write(PDF,'report.pdf'); z.write(HERE/'report.md','report.md')
        for f in (HERE/'figures').glob('*.png'): z.write(f,'figures/'+f.name)
        for f in ['evidence_manifest.json','build_report.py','make_figures.py','layout.json']: z.write(HERE/f,f)
        z.writestr('READ_ME.txt','Research evidence package for supervisory discussion, 8 November 2026.\nRead report.pdf first. report.md is editable.\nIncludes notebooks, application and research code snapshots, selected summaries and small raw inputs for report figures.\nThe complete raw NPZ datasets and many nested run files remain in the original repository.\nEmbedded experiment SHA manifests can reference files not included in this selected package.\nThe manifest uses repository-relative source paths; archive layout flattens the report and figures for convenience.\nNo new simulations were run for report preparation.\n')
    print('Evidence bundle',bundle,'files',len(files))


if __name__=='__main__':
    build(); evidence()
