#!/usr/bin/env python3
"""Render AlphaFold cartoons with 3Dmol.js in headless chromium.
All views share one camera frame (BCAT1's), since BCAT2 is superposed on BCAT1."""
import json, base64, os, sys
from playwright.sync_api import sync_playwright

D = "/root/bc/data"; R = "/root/bc/render"; OUT = "/root/bc/out/panels"
os.makedirs(OUT, exist_ok=True)
W, H = 1400, 1400          # css px; deviceScaleFactor 3 -> 4200 px PNG

BCAT1 = open(f"{D}/AF_P54687.pdb").read()
BCAT2 = open(f"{D}/AF_O15382_on_BCAT1.pdb").read()

GREEN, PURPLE, GREY = "#1b9e77", "#7570b3", "#c9ccd1"
# AlphaFold pLDDT palette
PL = [(90, "#0053D6"), (70, "#65CBF3"), (50, "#FFDB13"), (0, "#FF7D45")]

VIEWS = {
    # name: (js building models/styles)
    "bcat1_plddt": f"""
      var m = v.addModel(P1, "pdb");
      v.setStyle({{model:0}}, {{cartoon:{{colorfunc: plddt, thickness:0.4, arrows:true}}}});
    """,
    "bcat2_plddt": f"""
      var m0 = v.addModel(P1, "pdb"); v.setStyle({{model:0}}, {{}});
      var m1 = v.addModel(P2, "pdb");
      v.setStyle({{model:1}}, {{cartoon:{{colorfunc: plddt, thickness:0.4, arrows:true}}}});
    """,
    "superposition": f"""
      var m0 = v.addModel(P1, "pdb");
      var m1 = v.addModel(P2, "pdb");
      v.setStyle({{model:0}}, {{cartoon:{{color:"{GREY}", thickness:0.3, opacity:0.75, arrows:true}}}});
      v.setStyle({{model:1}}, {{cartoon:{{color:"{GREY}", thickness:0.3, opacity:0.75, arrows:true}}}});
      v.setStyle({{model:0, resi:"1-23"}}, {{cartoon:{{color:"{GREEN}", thickness:0.5, arrows:true}}}});
      v.setStyle({{model:1, resi:"1-30"}}, {{cartoon:{{color:"{PURPLE}", thickness:0.5, arrows:true}}}});
      v.addStyle({{model:0, resi:222}}, {{stick:{{color:"{GREEN}", radius:0.22}}}});
      v.addStyle({{model:1, resi:229}}, {{stick:{{color:"{PURPLE}", radius:0.22}}}});
      v.addStyle({{model:0, resi:[335,338]}}, {{stick:{{color:"#d95f02", radius:0.22}}}});
    """,
    "bcat1_nterm": f"""
      var m = v.addModel(P1, "pdb");
      v.setStyle({{model:0}}, {{cartoon:{{color:"{GREY}", thickness:0.35, arrows:true}}}});
      v.setStyle({{model:0, resi:"1-23"}}, {{cartoon:{{color:"{GREEN}", thickness:0.55, arrows:true}}}});
      v.addStyle({{model:0, resi:222}}, {{stick:{{color:"#222222", radius:0.22}}}});
      v.addStyle({{model:0, resi:[335,338]}}, {{stick:{{color:"#d95f02", radius:0.25}}}});
    """,
    "bcat2_nterm": f"""
      var m0 = v.addModel(P1, "pdb"); v.setStyle({{model:0}}, {{}});
      var m1 = v.addModel(P2, "pdb");
      v.setStyle({{model:1}}, {{cartoon:{{color:"{GREY}", thickness:0.35, arrows:true}}}});
      v.setStyle({{model:1, resi:"1-30"}}, {{cartoon:{{color:"{PURPLE}", thickness:0.55, arrows:true}}}});
      v.addStyle({{model:1, resi:229}}, {{stick:{{color:"#222222", radius:0.22}}}});
    """,
}

HTML = """<!doctype html><html><head><meta charset="utf-8">
<script src="3Dmol-min.js"></script>
<style>html,body{margin:0;background:#fff}#gl{width:%dpx;height:%dpx;position:relative}</style>
</head><body><div id="gl"></div><script>
var P1 = %s, P2 = %s;
function plddt(atom){var b=atom.b; return b>90?"#0053D6":b>70?"#65CBF3":b>50?"#FFDB13":"#FF7D45";}
var v = $3Dmol.createViewer("gl", {backgroundColor:"white", antialias:true});
%s
v.zoomTo({model:0});
v.rotate(%f,"y"); v.rotate(%f,"x"); v.rotate(%f,"z");
v.zoom(%f);
v.render();
window.__done = true;
</script></body></html>
"""

ROT = (float(sys.argv[1]) if len(sys.argv) > 1 else -25,
       float(sys.argv[2]) if len(sys.argv) > 2 else 15,
       float(sys.argv[3]) if len(sys.argv) > 3 else 0)
ZOOM = float(sys.argv[4]) if len(sys.argv) > 4 else 1.25

with sync_playwright() as p:
    b = p.chromium.launch(args=["--use-gl=swiftshader", "--enable-webgl",
                                "--ignore-gpu-blocklist"])
    pg = b.new_page(viewport={"width": W, "height": H}, device_scale_factor=3)
    for name, body in VIEWS.items():
        html = HTML % (W, H, json.dumps(BCAT1), json.dumps(BCAT2), body,
                       ROT[0], ROT[1], ROT[2], ZOOM)
        open(f"{R}/{name}.html", "w").write(html)
        pg.goto(f"file://{R}/{name}.html")
        pg.wait_for_function("window.__done === true", timeout=60000)
        pg.wait_for_timeout(1500)
        pg.locator("#gl").screenshot(path=f"{OUT}/3d_{name}.png")
        print("rendered", name)
    b.close()
