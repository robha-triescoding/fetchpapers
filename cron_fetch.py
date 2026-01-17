import json
import os
import sys
from Bio import Entrez
from google import genai
from datetime import datetime
import time

# Use environment variables for secrets (set by GitHub Actions)
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL")
ENTREZ_API_KEY = os.getenv("ENTREZ_API_KEY")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")

# Safety check for environment variables
missing_vars = []
if not ENTREZ_EMAIL: missing_vars.append("ENTREZ_EMAIL")
if not ENTREZ_API_KEY: missing_vars.append("ENTREZ_API_KEY")
if not GEMINI_API_KEY: missing_vars.append("GEMINI_API_KEY")

if missing_vars:
    print(f"ERROR: Missing environment variables: {', '.join(missing_vars)}")
    print("Please ensure these are set in GitHub Repository Secrets (Settings > Secrets and variables > Actions).")
    exit(1)

Entrez.email = ENTREZ_EMAIL
Entrez.api_key = ENTREZ_API_KEY
client = genai.Client(api_key=GEMINI_API_KEY)

def run_automated_fetch(set_id=None, day_of_week=None):
    # 1. Load Config
    try:
        with open("config.json", "r") as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config.json: {e}")
        return

    search_sets = config_data.get("search_sets", [])
    
    # Filter sets to run
    sets_to_run = []
    if set_id:
        sets_to_run = [s for s in search_sets if s['id'] == set_id]
        if not sets_to_run:
            print(f"Error: Search set with ID '{set_id}' not found.")
            return
    elif day_of_week:
        sets_to_run = [s for s in search_sets if s.get('schedule_day') == day_of_week]
        if not sets_to_run:
            print(f"No search sets scheduled for {day_of_week}.")
            return
    else:
        # Default to all if nothing specified (for backward compatibility or manual run)
        sets_to_run = search_sets

    for s_set in sets_to_run:
        print(f"\n--- Starting fetch for set: {s_set['name']} ({s_set['id']}) ---")
        print(f"Query: {s_set['query']}")
        
        # 2. Search PubMed
        try:
            handle = Entrez.esearch(
                db="pubmed",
                term=s_set['query'],
                reldate=s_set['days_back'],
                datetype="pdat",
                retmax=s_set['max_results']
            )
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]
        except Exception as e:
            print(f"Error searching PubMed for set {s_set['id']}: {e}")
            continue
        
        if not id_list:
            print(f"No new papers found for {s_set['name']}.")
            continue

        # 3. Fetch Details
        try:
            ids = ",".join(id_list)
            handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
        except Exception as e:
            print(f"Error fetching details for set {s_set['id']}: {e}")
            continue
        
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
                
                # 4. Analyze with Gemini
                prompt = f"Summarize this abstract in 3 bullet points highlighting the specific mechanism. Bold human trials.\n\nAbstract:\n{abstract}"
                response = client.models.generate_content(model='gemini-2.5-flash', contents=prompt)
                analysis = response.text
                
                papers.append({
                    'title': title,
                    'journal': journal,
                    'abstract': abstract,
                    'link': doi,
                    'analysis': analysis,
                    'fetched_at': datetime.now().strftime('%Y-%m-%d %H:%M')
                })
                print(f"Analyzed: {title[:50]}...")
                time.sleep(2) # Avoid hitting rate limits
            except Exception as e:
                print(f"Error processing paper: {e}")
                continue

        # 5. Save Results for this specific set
        results_filename = f"results_{s_set['id']}.json"
        with open(results_filename, "w") as f:
            json.dump(papers, f, indent=4)
        print(f"Successfully saved {len(papers)} papers to {results_filename}")

if __name__ == "__main__":
    # Handle arguments: python cron_fetch.py [set_id] [day_of_week]
    target_id = sys.argv[1] if len(sys.argv) > 1 else None
    target_day = sys.argv[2] if len(sys.argv) > 2 else None
    
    run_automated_fetch(set_id=target_id, day_of_week=target_day)
