"""
Messy Lab assistant in the browser.

    python3 assistant/app.py            # then open http://localhost:8765

Each question streams its database queries as they run, then the answer with the samples, files and
issue records behind it. Uses assistant.py for the work; the database is opened read-only.
"""
import json
import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assistant import ask  # noqa: E402

PAGE = (Path(__file__).resolve().parent / "page.html").read_text()
PORT = int(os.environ.get("PORT", 8765))


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def do_GET(self):
        u = urlparse(self.path)
        if u.path == "/":
            body = PAGE.encode()
            self.send_response(200); self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)
        elif u.path == "/ask":
            q = parse_qs(u.query).get("q", [""])[0].strip()
            self.send_response(200); self.send_header("Content-Type", "text/event-stream")
            self.send_header("Cache-Control", "no-cache"); self.end_headers()

            def send(kind, data):
                self.wfile.write(f"event: {kind}\ndata: {json.dumps(data, default=str)}\n\n".encode()); self.wfile.flush()
            try:
                ask(q, send) if q else send("error", {"error": "Type a question first."})
            except Exception as e:
                send("error", {"error": f"{type(e).__name__}: {e}"})
        else:
            self.send_response(404); self.end_headers()


if __name__ == "__main__":
    print(f"Messy Lab assistant on http://localhost:{PORT}")
    ThreadingHTTPServer(("127.0.0.1", PORT), Handler).serve_forever()
