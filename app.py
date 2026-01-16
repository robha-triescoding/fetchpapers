import streamlit as st
from Bio import Entrez
from datetime import datetime
from google import genai
import time
from fpdf import FPDF
import io

# --- STREAMLIT CONFIGURATION ---
st.set_page_config(page_title="PubMed Research Assistant", page_icon="🔬", layout="wide")

# --- SECRETS MANAGEMENT ---
try:
    ENTREZ_EMAIL = st.secrets["ENTREZ_EMAIL"]
    ENTREZ_API_KEY = st.secrets["ENTREZ_API_KEY"]
    GEMINI_API_KEY = st.secrets["GEMINI_API_KEY"]
    
    Entrez.email = ENTREZ_EMAIL
    Entrez.api_key = ENTREZ_API_KEY
    client = genai.Client(api_key=GEMINI_API_KEY)
except Exception as e:
    st.error("Missing configuration! Please ensure ENTREZ_EMAIL, ENTREZ_API_KEY, and GEMINI_API_KEY are set in Streamlit secrets.")
    st.stop()

# --- FUNCTIONS ---

def search_pubmed(query, days_back=7, max_results=5):
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            reldate=days_back,
            datetype="pdat",
            retmax=max_results
        )
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        st.error(f"Error searching PubMed: {e}")
        return []

def fetch_details(id_list):
    if not id_list:
        return []
    try:
        ids = ",".join(id_list)
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        papers = []
        for article in records['PubmedArticle']:
            try:
                article_data = article['MedlineCitation']['Article']
                title = article_data.get('ArticleTitle', 'No Title')
                journal = article_data.get('Journal', {}).get('Title', 'Unknown Journal')
                
                abstract_raw = article_data.get('Abstract', {}).get('AbstractText', [])
                abstract = " ".join(abstract_raw) if isinstance(abstract_raw, list) else str(abstract_raw)
                
                article_ids = article['PubmedData']['ArticleIdList']
                doi = next((f"https://doi.org/{item}" for item in article_ids if item.attributes.get('IdType') == 'doi'), "")
                
                papers.append({
                    'title': title,
                    'journal': journal,
                    'abstract': abstract,
                    'link': doi
                })
            except Exception:
                continue
        return papers
    except Exception as e:
        st.error(f"Error fetching details: {e}")
        return []

def analyze_abstract_with_retry(abstract, max_retries=3):
    if not abstract or abstract == "No Abstract Available":
        return "No abstract available to analyze."
    
    prompt = f"""
    You are an expert in metabolic disorders. 
    Summarize this abstract in 3 bullet points highlighting the specific mechanism discovered. 
    If the paper mentions human trials, bold that text.

    Abstract:
    {abstract}
    """
    
    for attempt in range(max_retries):
        try:
            # Using gemini-2.5-flash as requested
            response = client.models.generate_content(
                model='gemini-2.5-flash', 
                contents=prompt
            )
            return response.text
        except Exception as e:
            error_msg = str(e)
            if "429" in error_msg:
                if "PerDay" in error_msg:
                    return "DAILY QUOTA EXHAUSTED"
                time.sleep(30)
            else:
                return f"Error: {error_msg}"
    return "Failed to analyze after multiple retries."

def create_pdf(paper_details):
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    epw = pdf.epw
    
    pdf.set_font("helvetica", 'B', 16)
    date_str = datetime.now().strftime('%Y-%m-%d')
    pdf.cell(epw, 10, text=f"LATEST RESEARCH FINDINGS ({date_str})", align='C', new_x="LMARGIN", new_y="NEXT")
    pdf.ln(10)
    
    for i, paper in enumerate(paper_details, 1):
        pdf.set_x(pdf.l_margin)
        pdf.set_font("helvetica", 'B', 12)
        pdf.cell(epw, 10, text=f"PAPER #{i}", new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font("helvetica", 'B', 11)
        title_text = f"TITLE: {paper['title']}".encode('latin-1', 'replace').decode('latin-1')
        pdf.multi_cell(epw, 8, text=title_text, new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font("helvetica", 'I', 10)
        journal_text = f"JOURNAL: {paper.get('journal', 'Unknown Journal')}".encode('latin-1', 'replace').decode('latin-1')
        pdf.multi_cell(epw, 7, text=journal_text, new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font("helvetica", 'I', 9)
        pdf.set_text_color(0, 0, 255)
        link_text = f"LINK: {paper['link']}".encode('latin-1', 'replace').decode('latin-1')
        pdf.multi_cell(epw, 7, text=link_text, new_x="LMARGIN", new_y="NEXT")
        pdf.set_text_color(0, 0, 0)
        
        pdf.set_font("helvetica", 'B', 10)
        pdf.cell(epw, 8, text="GEMINI ANALYSIS:", new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font("helvetica", '', 10)
        analysis_text = paper.get('analysis', 'No analysis available.').replace('**', '')
        analysis_text = analysis_text.encode('latin-1', 'replace').decode('latin-1')
        pdf.multi_cell(epw, 6, text=analysis_text, new_x="LMARGIN", new_y="NEXT")
        
        pdf.ln(2)
        pdf.set_font("helvetica", 'B', 9)
        pdf.cell(epw, 7, text="ORIGINAL ABSTRACT:", new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font("helvetica", '', 8)
        pdf.set_text_color(50, 50, 50)
        orig_abstract = paper['abstract'].encode('latin-1', 'replace').decode('latin-1')
        pdf.multi_cell(epw, 5, text=orig_abstract, new_x="LMARGIN", new_y="NEXT")
        pdf.set_text_color(0, 0, 0)
        
        pdf.ln(5)
        pdf.line(pdf.l_margin, pdf.get_y(), pdf.l_margin + epw, pdf.get_y())
        pdf.ln(5)
    
    return pdf.output()

# --- UI LAYOUT ---

st.title("🔬 PubMed Research Assistant")
st.markdown("Fetch the latest research and get AI-powered summaries instantly.")

with st.sidebar:
    st.header("Search Parameters")
    search_query = st.text_area("Search Query", value='(("Adipose Tissue"[Title/Abstract] OR "Adipocytes"[Title/Abstract]) AND "Obesity"[Title/Abstract]) AND hasabstract[text]')
    days_back = st.slider("Days Back", 1, 30, 7)
    max_results = st.slider("Max Results", 1, 20, 5)
    fetch_button = st.button("Fetch & Analyze", type="primary")

if fetch_button:
    with st.spinner("Searching PubMed..."):
        ids = search_pubmed(search_query, days_back, max_results)
        
    if ids:
        with st.spinner(f"Fetching details for {len(ids)} papers..."):
            papers = fetch_details(ids)
            
        if papers:
            analyzed_papers = []
            for i, paper in enumerate(papers):
                with st.status(f"Analyzing Paper {i+1}/{len(papers)}...", expanded=True) as status:
                    st.write(f"**Title:** {paper['title']}")
                    analysis = analyze_abstract_with_retry(paper['abstract'])
                    paper['analysis'] = analysis
                    analyzed_papers.append(paper)
                    status.update(label=f"Analysis Complete for Paper {i+1}", state="complete")
            
            st.success("All papers analyzed!")
            
            # Download PDF Button
            pdf_data = create_pdf(analyzed_papers)
            st.download_button(
                label="📥 Download Results as PDF",
                data=pdf_data,
                file_name=f"fetchedpapers_{datetime.now().strftime('%Y%m%d')}.pdf",
                mime="application/pdf"
            )
            
            # Display Cards
            st.divider()
            for paper in analyzed_papers:
                with st.container(border=True):
                    st.subheader(paper['title'])
                    col1, col2 = st.columns([1, 4])
                    with col1:
                        st.caption(f"**Journal:**\n{paper['journal']}")
                        if paper['link']:
                            st.link_button("View Paper", paper['link'])
                    with col2:
                        st.markdown("**AI Summary:**")
                        st.write(paper['analysis'])
                    
                    with st.expander("Show Original Abstract"):
                        st.write(paper['abstract'])
    else:
        st.warning("No papers found matching your criteria.")
