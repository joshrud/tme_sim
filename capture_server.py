"""Serve the viewer and accept rendered frames back from the browser.

The 3D view only exists inside WebGL, so frames have to come from the browser. This serves
viewer/ as usual and adds POST /frame/<n>, which writes a PNG data URL to out/frames/.

  python capture_server.py [port] [outdir]
"""
import base64
import os
import sys
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "out/frames"
os.makedirs(OUTDIR, exist_ok=True)


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory="viewer", **k)

    def do_POST(self):
        if not self.path.startswith("/frame/"):
            self.send_error(404)
            return
        n = self.path.rsplit("/", 1)[1]
        raw = self.rfile.read(int(self.headers.get("Content-Length", 0))).decode()
        png = base64.b64decode(raw.split(",", 1)[1] if "," in raw else raw)
        with open(os.path.join(OUTDIR, f"f{int(n):04d}.png"), "wb") as fh:
            fh.write(png)
        self.send_response(204)
        self.end_headers()

    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-store")   # always serve the current viewer.js
        super().end_headers()

    def log_message(self, *a):
        pass


print(f"serving viewer/ on :{PORT}, frames -> {OUTDIR}/")
ThreadingHTTPServer(("127.0.0.1", PORT), Handler).serve_forever()
