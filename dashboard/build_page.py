"""
Messy Lab dashboard, step 2 of 2: put dashboard/dashboard_data.json into dashboard/template.html.
Writes dashboard/messy_lab_dashboard.html (the published page body) and dashboard/index.html (opens in any browser).

    python3 dashboard/build_data.py && python3 dashboard/build_page.py
"""
from pathlib import Path
HERE = Path(__file__).resolve().parent
data = (HERE / "dashboard_data.json").read_text().replace("</", "<\\/")
page = (HERE / "template.html").read_text().replace("/*__DATA__*/", data)
(HERE / "messy_lab_dashboard.html").write_text(page)
(HERE / "index.html").write_text('<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">\n'
                                 + page.replace("<style>", "<style>\n:root{color-scheme:light}", 1) + "\n</html>\n")
print(f"messy_lab_dashboard.html and index.html written ({len(page)/1e6:.1f} MB)")
