from bs4 import BeautifulSoup
import argparse
import os

def convert_rna_to_fixed_text(html_content, output_file):
    soup = BeautifulSoup(html_content, 'html.parser')
    blacklist = ["RNA STRUCTURE ANALYSIS", "STRUCTURAL INFORMATION"]
    std_w = 18
    seq_w = 45  # Width for Sequence column
    str_w = 60  # Width for Structure column

    with open(output_file, 'w', encoding='utf-8') as f:
        # 1. Header Information
        info_tip = soup.find('div', class_='info-tip')
        if info_tip:
            f.write("="*120 + "\n")
            f.write(" RNA ANALYZER RESULTS\n")
            f.write("="*120 + "\n")
            f.write(info_tip.get_text(separator=' | ', strip=True) + "\n\n")

        # 2. BASIC INFORMATION (first table in main)
        main_content = soup.find('main')
        first_table = main_content.find('table', class_='table-result') if main_content else None
        
        if first_table:
            f.write("[ BASIC INFORMATION ]\n")
            f.write("-" * 21 + "\n")
            for row in first_table.find_all('tr'):
                cols = row.find_all(['th', 'td'])
                if len(cols) >= 2:
                    f.write(f"{cols[0].get_text(strip=True).ljust(25)}: {cols[1].get_text(strip=True)}\n")
            f.write("\n")

        # 3. Process Analysis Boxes
        for box in soup.find_all('div', class_='box'):
            header_div = box.find('div', class_='box-header')
            if not header_div:
                continue
            
            section_name = header_div.get_text(strip=True).upper()
            
            # Check blacklist
            if any(b in section_name for b in blacklist):
                # print(f"Skipping blacklisted section: {section_name}")
                continue
            
            f.write(f"[ {section_name} ]\n")
            f.write("-" * (len(section_name) + 4) + "\n")

            # Get the box content div
            box_content = box.find('div', class_=['box-content', 'box-content-structure'])
            
            if box_content:
                # Capture standalone text/elements (not in tables or info divs)
                for element in box_content.children:
                    # Skip tables and info divs (we process these separately)
                    if element.name == 'table' or (element.name == 'div' and 
                        any(c in element.get('class', []) for c in ['info-warning', 'info-info', 'box'])):
                        continue
                    
                    # Get text from direct text nodes and simple elements
                    if element.name is None:  # Text node
                        text = str(element).strip()
                        if text:
                            f.write(f"{text}\n")
                    elif element.name in ['i', 'b', 'strong', 'em', 'p', 'br']:
                        if element.name == 'br':
                            continue  # Skip br tags
                        text = element.get_text(strip=True)
                        if text:
                            f.write(f"{text}\n")
            
            # Write warnings/info messages
            has_notes = False
            for note in box.find_all('div', class_=['info-warning', 'info-info']):
                if note.find_parent('div', class_='box') == box:
                    f.write(f"! {note.get_text(separator=' ', strip=True)}\n")
                    has_notes = True
            
            if has_notes:
                f.write("\n")  # Add line gap after info tags

            # Process tables
            for table in box.find_all('table', class_='table-result'):
                # Skip if this is the first table we already processed
                if table == first_table:
                    continue
                
                # Skip if table belongs to a nested box
                if table.find_parent('div', class_='box') != box:
                    continue

                rows = table.find_all('tr')
                if not rows:
                    continue

                # Get headers and filter out 'View'
                header_cells = rows[0].find_all(['th', 'td'])
                raw_headers = [h.get_text(strip=True) for h in header_cells]
                
                # Filter indices where we keep the column
                keep_indices = [i for i, h in enumerate(raw_headers) 
                               if "View" not in h and h != ""]
                headers = [raw_headers[i] for i in keep_indices]

                if len(headers) > 1:
                    # Multi-column table
                    # print(f"Processing table with headers: {headers}")
                    
                    # Print headers
                    header_line = ""
                    for h in headers:
                        if "Structure" in h:
                            w = str_w
                        elif "Sequence" in h:
                            w = seq_w
                        else:
                            w = std_w
                        header_line += h.ljust(w)
                    f.write(header_line + "\n")
                    f.write("." * len(header_line) + "\n")

                    # Build grid for rowspan handling
                    num_rows = len(rows)
                    num_cols = len(raw_headers)
                    grid = [["" for _ in range(num_cols)] for _ in range(num_rows)]

                    # Fill grid
                    for r_idx, row in enumerate(rows):
                        c_idx = 0
                        cells = row.find_all(['td', 'th'], recursive=False)
                        
                        for cell in cells:
                            # Find next empty column
                            while c_idx < num_cols and grid[r_idx][c_idx] != "":
                                c_idx += 1
                            
                            if c_idx >= num_cols:
                                break
                            
                            val = cell.get_text(strip=True)
                            rs = int(cell.get('rowspan', 1))
                            cs = int(cell.get('colspan', 1))
                            
                            # Fill grid cells
                            for r_offset in range(rs):
                                for c_offset in range(cs):
                                    nr = r_idx + r_offset
                                    nc = c_idx + c_offset
                                    if nr < num_rows and nc < num_cols:
                                        if r_offset == 0 and c_offset == 0:
                                            grid[nr][nc] = val
                                        else:
                                            grid[nr][nc] = " "
                            
                            c_idx += cs

                    # Print data rows
                    for r_idx in range(1, num_rows):
                        line = ""
                        for idx, c_idx in enumerate(keep_indices):
                            h_name = headers[idx]
                            if "Structure" in h_name:
                                w = str_w
                            elif "Sequence" in h_name:
                                w = seq_w
                            else:
                                w = std_w
                            line += grid[r_idx][c_idx].ljust(w)
                        f.write(line.rstrip() + "\n")
                
                else:
                    # Key-value table
                    # print(f"Processing key-value table")
                    for r in rows:
                        cols = r.find_all(['td', 'th'])
                        if len(cols) >= 2:
                            k = cols[0].get_text(strip=True)
                            v = cols[1].get_text(strip=True)
                            if "View" not in v and k:
                                f.write(f"{k.ljust(25)}: {v}\n")
            
            f.write("\n")

# Run the conversion
def main():
    # 1. Initialize the parser
    parser = argparse.ArgumentParser(description="Convert RNA HTML results to fixed text.")

    # 2. Add arguments
    # We can pass the job ID to construct the path, or the full path directly
    parser.add_argument('--input', '-i', type=str, required=True, 
                        help="Path to the input HTML file")
    
    parser.add_argument('--output', '-o', type=str, default='rna_final_report.txt', 
                        help="Path for the output text file (default: rna_final_report.txt)")

    # 3. Parse the arguments
    args = parser.parse_args()

    # 4. Use the variables in your logic
    if not os.path.exists(args.input):
        print(f"Error: The file {args.input} does not exist.")
        return

    try:
        with open(args.input, 'r', encoding='utf-8') as f:
            html_content = f.read()
        
        convert_rna_to_fixed_text(html_content, args.output)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
    
    
