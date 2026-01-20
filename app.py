import streamlit as st
from Bio import Entrez
from datetime import datetime
from google import genai
import time
from fpdf import FPDF
import io
import json
import os
import requests
import base64

# --- STREAMLIT CONFIGURATION ---
st.set_page_config(page_title="PubMed Research Assistant", page_icon="🔬", layout="wide")

# --- SECRETS & CONFIG ---
try:
    ENTREZ_EMAIL = st.secrets["ENTREZ_EMAIL"]
    ENTREZ_API_KEY = st.secrets["ENTREZ_API_KEY"]
    GEMINI_API_KEY = st.secrets["GEMINI_API_KEY"]
    # Admin Password for sharing
    ADMIN_PASSWORD = st.secrets.get("ADMIN_PASSWORD", "research2026")
    # Optional for saving settings back to GitHub
    GITHUB_TOKEN = st.secrets.get("GITHUB_TOKEN")
    GITHUB_REPO = st.secrets.get("GITHUB_REPO") # format: "username/repo"
    
    Entrez.email = ENTREZ_EMAIL
    Entrez.api_key = ENTREZ_API_KEY
    client = genai.Client(api_key=GEMINI_API_KEY)
except Exception as e:
    st.error("Missing configuration! Please ensure ENTREZ_EMAIL, ENTREZ_API_KEY, and GEMINI_API_KEY are set in Streamlit secrets.")
    st.stop()

# Load default config from file
def load_config():
    try:
        with open("config.json", "r", encoding='utf-8') as f:
            config = json.load(f)
            # Migration for old config format
            if "search_query" in config:
                return {
                    "search_sets": [
                        {
                            "id": "default",
                            "name": "Default Set",
                            "query": config["search_query"],
                            "days_back": config.get("days_back", 7),
                            "max_results": config.get("max_results", 5),
                            "schedule_day": "Sunday"
                        }
                    ]
                }
            return config
    except:
        return {
            "search_sets": [
                {
                    "id": "adipose_obesity",
                    "name": "Adipose Tissue in Obesity",
                    "query": "((\"Adipose Tissue\"[Title/Abstract] OR \"Adipocytes\"[Title/Abstract]) AND \"Obesity\"[Title/Abstract]) AND hasabstract[text]",
                    "days_back": 7,
                    "max_results": 5,
                    "schedule_day": "Sunday"
                }
            ]
        }

def build_pubmed_query(terms):
    """
    Constructs a PubMed query string from a list of terms and operators.
    Each term is a dict: {"text": "...", "operator": "AND/OR/NOT", "field": "Title/Abstract/None"}
    """
    if not terms:
        return ""
    
    query = ""
    for i, term in enumerate(terms):
        text = term['text'].strip()
        if not text: continue
        
        # Add quotes if multi-word
        if " " in text and not text.startswith('"'):
            text = f'"{text}"'
            
        field_suffix = f"[{term['field']}]" if term['field'] != "None" else ""
        
        if i == 0:
            query = f"{text}{field_suffix}"
        else:
            op = term.get('operator', 'AND')
            query = f"({query} {op} {text}{field_suffix})"
            
    # Always ensure it has an abstract
    if "hasabstract[text]" not in query:
        query = f"({query}) AND hasabstract[text]"
        
    return query

# Save config to GitHub (since Streamlit Cloud is read-only)
def save_config_to_github(new_config):
    if not GITHUB_TOKEN or not GITHUB_REPO:
        st.warning("GITHUB_TOKEN and GITHUB_REPO not found in secrets. Settings saved for this session only.")
        return False
    
    url = f"https://api.github.com/repos/{GITHUB_REPO}/contents/config.json"
    headers = {"Authorization": f"token {GITHUB_TOKEN}", "Accept": "application/vnd.github.v3+json"}
    
    # Get current file sha
    r = requests.get(url, headers=headers)
    sha = r.json().get("sha") if r.status_code == 200 else None
    
    content = base64.b64encode(json.dumps(new_config, indent=4, ensure_ascii=False).encode('utf-8')).decode()
    data = {
        "message": "Update search settings via Streamlit UI",
        "content": content,
        "sha": sha
    }
    
    r = requests.put(url, headers=headers, json=data)
    return r.status_code in [200, 201]

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
    
    # Add Unicode-capable font
    font_name = "helvetica"
    try:
        # Use a more robust way to load fonts by downloading them to a local temp path
        # this avoids issues with fpdf2 not handling certain URL types or redirects
        fonts = {
            "Roboto-Regular.ttf": "https://raw.githubusercontent.com/google/fonts/main/apache/roboto/static/Roboto-Regular.ttf",
            "Roboto-Bold.ttf": "https://raw.githubusercontent.com/google/fonts/main/apache/roboto/static/Roboto-Bold.ttf",
            "Roboto-Italic.ttf": "https://raw.githubusercontent.com/google/fonts/main/apache/roboto/static/Roboto-Italic.ttf"
        }
        
        for name, url in fonts.items():
            if not os.path.exists(name):
                resp = requests.get(url, timeout=10)
                if resp.status_code == 200:
                    with open(name, "wb") as f:
                        f.write(resp.content)
        
        pdf.add_font("Roboto", style="", fname="Roboto-Regular.ttf")
        pdf.add_font("Roboto", style="B", fname="Roboto-Bold.ttf")
        pdf.add_font("Roboto", style="I", fname="Roboto-Italic.ttf")
        font_name = "Roboto"
    except Exception as e:
        # Fallback to Helvetica if font loading fails
        st.warning(f"Could not load Unicode font, falling back to Helvetica. Special characters may not display correctly. Error: {e}")
        font_name = "helvetica"

    pdf.add_page()
    epw = pdf.epw
    
    pdf.set_font(font_name, 'B', 16)
    date_str = datetime.now().strftime('%Y-%m-%d')
    pdf.cell(epw, 10, text=f"LATEST RESEARCH FINDINGS ({date_str})", align='C', new_x="LMARGIN", new_y="NEXT")
    pdf.ln(10)
    
    for i, paper in enumerate(paper_details, 1):
        # Function to clean text for the current font
        def clean_for_pdf(txt):
            if font_name == "helvetica":
                return txt.encode('latin-1', 'replace').decode('latin-1')
            return str(txt)

        pdf.set_x(pdf.l_margin)
        pdf.set_font(font_name, 'B', 12)
        pdf.cell(epw, 10, text=clean_for_pdf(f"PAPER #{i}"), new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font(font_name, 'B', 11)
        pdf.multi_cell(epw, 8, text=clean_for_pdf(f"TITLE: {paper['title']}"), new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font(font_name, 'I', 10)
        pdf.multi_cell(epw, 7, text=clean_for_pdf(f"JOURNAL: {paper.get('journal', 'Unknown Journal')}"), new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font(font_name, 'I', 9)
        pdf.set_text_color(0, 0, 255)
        pdf.multi_cell(epw, 7, text=clean_for_pdf(f"LINK: {paper['link']}"), new_x="LMARGIN", new_y="NEXT")
        pdf.set_text_color(0, 0, 0)
        
        pdf.set_font(font_name, 'B', 10)
        pdf.cell(epw, 8, text=clean_for_pdf("GEMINI ANALYSIS:"), new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font(font_name, '', 10)
        analysis_text = paper.get('analysis', 'No analysis available.').replace('**', '')
        pdf.multi_cell(epw, 6, text=clean_for_pdf(analysis_text), new_x="LMARGIN", new_y="NEXT")
        
        pdf.ln(2)
        pdf.set_font(font_name, 'B', 9)
        pdf.cell(epw, 7, text=clean_for_pdf("ORIGINAL ABSTRACT:"), new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font(font_name, '', 8)
        pdf.set_text_color(50, 50, 50)
        pdf.multi_cell(epw, 5, text=clean_for_pdf(paper['abstract']), new_x="LMARGIN", new_y="NEXT")
        pdf.set_text_color(0, 0, 0)
        
        pdf.ln(5)
        pdf.line(pdf.l_margin, pdf.get_y(), pdf.l_margin + epw, pdf.get_y())
        pdf.ln(5)
    
    return pdf.output()

# --- UI LAYOUT ---

st.title("🔬 PubMed Research Assistant")
st.markdown("Fetch the latest research and get AI-powered summaries instantly.")

# Load initial config
if 'config' not in st.session_state:
    st.session_state.config = load_config()

# Initialize session state for papers
if 'analyzed_papers' not in st.session_state:
    st.session_state.analyzed_papers = {} # Store as {set_id: [papers]}

# Sidebar for Admin and Configuration
with st.sidebar:
    st.header("Admin Controls")
    password_input = st.text_input("Enter Admin Password to modify settings", type="password")
    is_admin = (password_input == ADMIN_PASSWORD)
    fetch_button = False
    save_button = False
    
    if is_admin:
        st.success("Admin Access Granted")
        st.divider()
        
        st.header("Search Set Management")
        
        # Select or Create Set
        set_names = [s['name'] for s in st.session_state.config['search_sets']]
        selected_set_name = st.selectbox("Select Search Set to Edit", ["+ Add New Set"] + set_names)
        
        st.divider() # Added divider for clarity
        
        if selected_set_name == "+ Add New Set":
            new_set_name = st.text_input("New Set Name", placeholder="e.g. Redox Regulation")
            if st.button("Create Set", use_container_width=True):
                if new_set_name:
                    # Sanitize ID to be filename-safe (no slashes, spaces, or special characters)
                    import re
                    new_id = re.sub(r'[^a-zA-Z0-9_-]', '_', new_set_name.lower().replace(" ", "_"))
                    st.session_state.config['search_sets'].append({
                        "id": new_id,
                        "name": new_set_name,
                        "query": "",
                        "days_back": 7,
                        "max_results": 5,
                        "schedule_day": "Sunday"
                    })
                    # Save to GitHub immediately so it persists on reload
                    if save_config_to_github(st.session_state.config):
                        st.success(f"Set '{new_set_name}' created and saved!")
                        time.sleep(1)
                        st.rerun()
                    else:
                        st.error("Set created in session, but failed to save to GitHub.")
        else:
            # Edit existing set
            s_set = next(s for s in st.session_state.config['search_sets'] if s['name'] == selected_set_name)
            
            s_set['name'] = st.text_input("Set Name", value=s_set['name'])
            s_set['schedule_day'] = st.selectbox("Schedule Day", 
                ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"],
                index=["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"].index(s_set.get('schedule_day', 'Sunday')))
            
            st.divider()
            st.subheader("Query Builder")
            
            # Simple Query Builder State
            if f"builder_{s_set['id']}" not in st.session_state:
                st.session_state[f"builder_{s_set['id']}"] = []
                # Try to parse existing query if builder is empty? 
                # (Complex, let's just allow manual editing or builder starting from scratch)
            
            builder_terms = st.session_state[f"builder_{s_set['id']}"]
            
            # Show current terms
            for idx, term in enumerate(builder_terms):
                col_t, col_o, col_f, col_d = st.columns([3, 2, 2, 1])
                with col_t: st.text(term['text'])
                with col_o: st.text(term['operator'] if idx > 0 else "-")
                with col_f: st.text(term['field'])
                with col_d: 
                    if st.button("🗑️", key=f"del_{s_set['id']}_{idx}"):
                        builder_terms.pop(idx)
                        st.rerun()
            
            # Add new term
            with st.container(border=True):
                new_text = st.text_input("Keyword", key=f"key_{s_set['id']}")
                col_op, col_fi = st.columns(2)
                with col_op:
                    new_op = st.selectbox("Operator", ["AND", "OR", "NOT"], key=f"op_{s_set['id']}")
                with col_fi:
                    new_fi = st.selectbox("Field", ["Title/Abstract", "Title", "Abstract", "None"], key=f"fi_{s_set['id']}")
                
                if st.button("Add Term", use_container_width=True):
                    if new_text:
                        builder_terms.append({"text": new_text, "operator": new_op, "field": new_fi})
                        s_set['query'] = build_pubmed_query(builder_terms)
                        # Clear text input after adding
                        st.rerun()
            
            s_set['query'] = st.text_area("Final PubMed Query", value=s_set['query'], help="You can also manually edit this.")
            
            col1, col2 = st.columns(2)
            with col1:
                s_set['days_back'] = st.number_input("Days Back", 1, 365, s_set['days_back'])
            with col2:
                s_set['max_results'] = st.number_input("Max Results", 1, 50, s_set['max_results'])
            
            st.divider()
            col_fetch, col_save = st.columns(2)
            with col_fetch:
                fetch_button = st.button("Fetch Now", type="primary", use_container_width=True)
            with col_save:
                save_button = st.button("Save Config", use_container_width=True)
                
            if st.button("Delete Set", type="secondary", use_container_width=True):
                st.session_state.config['search_sets'] = [s for s in st.session_state.config['search_sets'] if s['id'] != s_set['id']]
                if save_config_to_github(st.session_state.config):
                    st.success("Set deleted and config updated!")
                    time.sleep(1)
                    st.rerun()
                else:
                    st.error("Deleted from session, but failed to update GitHub.")
    else:
        if password_input:
            st.error("Incorrect Password")
        st.info("The app is currently in **Read-Only Mode**. Select a research topic below to view the latest results.")

# --- MAIN AREA ---

# If admin saved
if is_admin and save_button:
    if save_config_to_github(st.session_state.config):
        st.sidebar.success("All settings saved to GitHub!")
    else:
        st.sidebar.error("Could not save to GitHub.")

# If admin clicked fetch
if is_admin and fetch_button:
    # Use the selected set from the edit menu
    s_set = next(s for s in st.session_state.config['search_sets'] if s['name'] == selected_set_name)
    
    with st.spinner(f"Searching PubMed for '{s_set['name']}'..."):
        ids = search_pubmed(s_set['query'], s_set['days_back'], s_set['max_results'])
        
    if ids:
        with st.spinner(f"Fetching details for {len(ids)} papers..."):
            papers = fetch_details(ids)
            
        if papers:
            st.session_state.analyzed_papers[s_set['id']] = [] # Reset this set
            for i, paper in enumerate(papers):
                with st.status(f"Analyzing Paper {i+1}/{len(papers)}...", expanded=False) as status:
                    st.write(f"**Title:** {paper['title']}")
                    analysis = analyze_abstract_with_retry(paper['abstract'])
                    paper['analysis'] = analysis
                    st.session_state.analyzed_papers[s_set['id']].append(paper)
                    status.update(label=f"Analysis Complete for Paper {i+1}", state="complete")
            
            st.success(f"Analysis for '{s_set['name']}' complete!")
    else:
        st.warning("No papers found matching your criteria.")

# Display Results Selection
available_sets = st.session_state.config['search_sets']
if available_sets:
    tabs = st.tabs([s['name'] for s in available_sets])
    
    for i, s_set in enumerate(available_sets):
        with tabs[i]:
            # Try to load papers for this set if not in session state
            if s_set['id'] not in st.session_state.analyzed_papers:
                # Sanitize ID to match the filename logic in cron_fetch.py
                safe_id = "".join([c if c.isalnum() or c in ('-', '_') else '_' for c in s_set['id']])
                try:
                    with open(f"results_{safe_id}.json", "r", encoding='utf-8') as f:
                        st.session_state.analyzed_papers[s_set['id']] = json.load(f)
                except:
                    # Fallback to the exact set ID if the safe_id doesn't exist
                    try:
                        with open(f"results_{s_set['id']}.json", "r", encoding='utf-8') as f:
                            st.session_state.analyzed_papers[s_set['id']] = json.load(f)
                    except:
                        # Final fallback for migration
                        if i == 0:
                            try:
                                with open("results.json", "r", encoding='utf-8') as f:
                                    st.session_state.analyzed_papers[s_set['id']] = json.load(f)
                            except:
                                st.session_state.analyzed_papers[s_set['id']] = []
                        else:
                            st.session_state.analyzed_papers[s_set['id']] = []
            
            papers = st.session_state.analyzed_papers[s_set['id']]
            
            if papers:
                col_title, col_dl = st.columns([3, 1])
                with col_title:
                    st.subheader(f"Latest Results for {s_set['name']}")
                with col_dl:
                    pdf_output = create_pdf(papers)
                    pdf_bytes = bytes(pdf_output)
                    st.download_button(
                        label="📥 Download PDF",
                        data=pdf_bytes,
                        file_name=f"research_{s_set['id']}_{datetime.now().strftime('%Y%m%d')}.pdf",
                        mime="application/pdf",
                        key=f"dl_{s_set['id']}"
                    )
                
                for paper in papers:
                    with st.container(border=True):
                        st.markdown(f"### {paper['title']}")
                        col1, col2 = st.columns([1, 4])
                        with col1:
                            st.caption(f"**Journal:**\n{paper['journal']}")
                            if paper['fetched_at']:
                                st.caption(f"**Fetched:**\n{paper['fetched_at']}")
                            if paper['link']:
                                st.link_button("View Paper", paper['link'])
                        with col2:
                            st.markdown("**AI Summary:**")
                            st.write(paper['analysis'])
                        
                        with st.expander("Show Original Abstract"):
                            st.write(paper['abstract'])
            else:
                st.info(f"No results found for {s_set['name']} yet. Automated fetch runs on {s_set.get('schedule_day', 'Sunday')}s.")
else:
    st.info("No search sets configured. Please enter Admin mode to add one.")
