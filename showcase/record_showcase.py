"""Record a captioned showcase video of Protein Analyzer (for LinkedIn etc.).

Needs: pip install playwright imageio-ffmpeg && python -m playwright install chromium
Run the app first (python app.py), then:  python showcase/record_showcase.py
Output: showcase/protein_analyzer_showcase.mp4 (1280x720) and showcase/screenshots/*.png
"""
from __future__ import annotations

import asyncio
import json
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path

import imageio_ffmpeg
from playwright.async_api import async_playwright

BASE = "http://127.0.0.1:5000"
OUT = Path(__file__).resolve().parent
CLIPS, SHOTS = OUT / "clips", OUT / "screenshots"
W, H = 1280, 720

# Visible cursor + click ripple + caption bar, injected into every page.
INIT_JS = r"""
(() => {
  const add = () => {
    if (document.getElementById('fakecursor')) return;
    const c = document.createElement('div'); c.id = 'fakecursor';
    c.style.cssText = 'position:fixed;left:-50px;top:-50px;width:22px;height:22px;z-index:99999;pointer-events:none;transition:transform .08s';
    c.innerHTML = '<svg width="22" height="22" viewBox="0 0 22 22"><path d="M3 2 L3 18 L7.5 13.5 L10.5 20 L13 19 L10 12.5 L16 12.5 Z" fill="#111" stroke="#fff" stroke-width="1.5"/></svg>';
    document.body.appendChild(c);
    addEventListener('mousemove', e => { c.style.left = e.clientX + 'px'; c.style.top = e.clientY + 'px'; }, true);
    addEventListener('mousedown', e => {
      const r = document.createElement('div');
      r.style.cssText = `position:fixed;left:${e.clientX - 18}px;top:${e.clientY - 18}px;width:36px;height:36px;border-radius:50%;border:3px solid #2f5bea;z-index:99998;pointer-events:none;opacity:.9;transition:all .45s ease-out`;
      document.body.appendChild(r); requestAnimationFrame(() => { r.style.transform = 'scale(1.8)'; r.style.opacity = '0'; });
      setTimeout(() => r.remove(), 500);
    }, true);
  };
  if (document.body) add(); else addEventListener('DOMContentLoaded', add);
  window.__caption = (t, sub) => {
    let b = document.getElementById('capbar');
    if (!b) { b = document.createElement('div'); b.id = 'capbar';
      b.style.cssText = 'position:fixed;left:50%;bottom:26px;transform:translateX(-50%);max-width:88%;background:rgba(15,23,42,.88);color:#fff;padding:12px 22px;border-radius:14px;font:600 24px/1.3 system-ui,Segoe UI,sans-serif;z-index:99997;box-shadow:0 8px 30px rgba(0,0,0,.3);text-align:center;transition:opacity .3s;pointer-events:none';
      document.body.appendChild(b); }
    b.innerHTML = t + (sub ? `<div style="font:400 17px/1.35 system-ui,Segoe UI,sans-serif;opacity:.85;margin-top:3px">${sub}</div>` : '');
    b.style.opacity = t ? '1' : '0';
  };
})();
"""

CARD_HTML = """<!doctype html><html><body style="margin:0;height:100vh;display:grid;place-items:center;
background:linear-gradient(135deg,#1b2f7a,#2f5bea 55%,#0ea5b7);font-family:system-ui,Segoe UI,sans-serif;color:#fff">
<div style="text-align:center;max-width:1000px;padding:0 40px">{body}</div></body></html>"""


def warm(paths_get=(), posts=()):
    """Call the API once before recording so scenes don't show loading spinners."""
    for p in paths_get:
        try:
            urllib.request.urlopen(BASE + p, timeout=300).read()
        except Exception as e:
            print("warm failed", p, e)
    for p, body in posts:
        try:
            req = urllib.request.Request(BASE + p, json.dumps(body).encode(), {"Content-Type": "application/json"})
            urllib.request.urlopen(req, timeout=300).read()
        except Exception as e:
            print("warm failed", p, e)


class Scene:
    def __init__(self, browser, name):
        self.browser, self.name, self.cuts = browser, name, []

    async def __aenter__(self):
        self.ctx = await self.browser.new_context(viewport={"width": W, "height": H}, color_scheme="light",
                                                  record_video_dir=str(CLIPS / self.name),
                                                  record_video_size={"width": W, "height": H})
        await self.ctx.add_init_script(INIT_JS)
        self.page = await self.ctx.new_page()
        self.t0 = time.monotonic()
        return self

    async def __aexit__(self, *exc):
        await self.ctx.close()
        vids = list((CLIPS / self.name).glob("*.webm"))
        if vids:
            vids[0].replace(CLIPS / f"{self.name}.webm")
            shutil.rmtree(CLIPS / self.name, ignore_errors=True)
        (CLIPS / f"{self.name}.cuts.json").write_text(json.dumps(self.cuts))

    def now(self):
        return time.monotonic() - self.t0

    async def cap(self, text, sub=None):
        await self.page.evaluate("([t, s]) => window.__caption && window.__caption(t, s)", [text, sub])

    async def wait(self, ms):
        await self.page.wait_for_timeout(ms)

    async def move_to(self, selector=None, x=None, y=None, steps=25):
        if selector:
            box = await self.page.locator(selector).first.bounding_box()
            x, y = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
        await self.page.mouse.move(x, y, steps=steps)
        return x, y

    async def click(self, selector, pause=350):
        await self.move_to(selector)
        await self.wait(pause)
        await self.page.locator(selector).first.click()

    async def scroll_to(self, selector, block="center", wait=900):
        await self.page.evaluate("([s, b]) => document.querySelector(s).scrollIntoView({behavior: 'smooth', block: b})", [selector, block])
        await self.wait(wait)

    async def drag(self, selector, dx, dy, steps=40):
        box = await self.page.locator(selector).first.bounding_box()
        x, y = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
        await self.page.mouse.move(x, y, steps=10)
        await self.page.mouse.down()
        await self.page.mouse.move(x + dx, y + dy, steps=steps)
        await self.page.mouse.up()

    async def shot(self, name):
        SHOTS.mkdir(exist_ok=True)
        await self.page.evaluate("() => { const b = document.getElementById('capbar'); if (b) b.style.visibility = 'hidden'; const c = document.getElementById('fakecursor'); if (c) c.style.visibility = 'hidden'; }")
        await self.page.screenshot(path=str(SHOTS / f"{name}.png"))
        await self.page.evaluate("() => { const b = document.getElementById('capbar'); if (b) b.style.visibility = ''; const c = document.getElementById('fakecursor'); if (c) c.style.visibility = ''; }")


async def card(browser, name, body, seconds):
    async with Scene(browser, name) as s:
        await s.page.set_content(CARD_HTML.format(body=body))
        await s.wait(int(seconds * 1000))


async def record():
    CLIPS.mkdir(parents=True, exist_ok=True)
    warm(paths_get=["/api/search?q=hemoglobin%20AND%20organism_id%3A9606%20AND%20reviewed%3Atrue&size=25",
                    "/api/protein/P68871", "/api/structure/P68871", "/api/protein/P68871/structure-features",
                    "/api/protein/P68871/interactions?limit=15", "/api/protein/P68871/pathways",
                    "/api/protein/P68871/domains", "/api/protein/P68871/alphamissense", "/api/protein/P68871/variants"],
         posts=[("/api/mutate", {"accession": "P68871", "mutation": "E7V"}),
                ("/api/msa", {"text": "P69905 P68871 P02042 P02100 P69891 P02144"}),
                ("/api/superpose", {"a": "P69905", "b": "P02144"}),
                ("/api/classify", {"text": "P00533", "backend": "esm2"})])
    async with async_playwright() as p:
        browser = await p.chromium.launch(args=["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"])

        await card(browser, "00_title", """<div style="font-size:72px;margin-bottom:6px">🧬</div>
          <div style="font-size:60px;font-weight:800;letter-spacing:-1px">Protein Analyzer</div>
          <div style="font-size:26px;opacity:.92;margin-top:14px">Explore any protein with real scientific data:<br>
          3D structures · disease mutations · evolution · machine learning</div>""", 3.5)

        # 0. Home: new design, dark mode, quick-jump palette
        async with Scene(browser, "00b_home") as s:
            pg = s.page
            await pg.goto(BASE + "/")
            await s.cap("One friendly home for every tool", "search by name, accession or sequence")
            await s.move_to("#homeq"); await s.wait(1200)
            await s.move_to(".tools a[href='/mutate']"); await s.wait(900)
            await s.move_to(".tools a[href='/classify']"); await s.wait(900)
            await s.shot("00_home")
            await s.cap("Light or dark, your choice", "or follow your computer's setting")
            await s.click("[data-theme-set='dark']"); await s.wait(2200)
            await s.shot("00_home_dark")
            await s.cap("Press / anywhere to jump to a protein or tool")
            await pg.keyboard.press("/"); await s.wait(500)
            await pg.keyboard.type("mutation", delay=110); await s.wait(1600)
            await pg.keyboard.press("Escape"); await s.wait(300)
            await s.click("[data-theme-set='light']"); await s.wait(700)

        # 1. Search -> protein page
        async with Scene(browser, "01_search") as s:
            pg = s.page
            await pg.goto(BASE + "/search")
            await s.cap("Search UniProt with simple filters", "250 million protein records, no query syntax needed")
            await s.click("#q"); await pg.keyboard.type("hemoglobin", delay=90)
            await s.move_to("#org"); await pg.select_option("#org", "organism_id:9606")
            await pg.evaluate("document.querySelector('#org').dispatchEvent(new Event('input', {bubbles: true}))")
            await s.wait(500)
            await s.click("#go")
            await pg.wait_for_selector("#rows tr", timeout=60000)
            await s.wait(1500)
            await s.shot("01_search")
            await s.click("#rows tr[data-acc='P68871'] a")
            await pg.wait_for_selector("#out h1", timeout=60000)
            await s.cap("Every protein in one place", "function, gene ontology, enzyme class, keywords…")
            await s.wait(2200)

        # 2. Protein page: 3D, secondary structure, network
        async with Scene(browser, "02_protein") as s:
            pg = s.page
            await pg.goto(BASE + "/protein/P68871")
            await pg.wait_for_selector("#viewer canvas", timeout=60000)
            await pg.wait_for_selector("#sstrack svg", timeout=60000)
            await s.cap("Save anything to a project", "results, notes and data sources kept together")
            await s.click("#out .save-btn")
            await pg.wait_for_selector("#savedlg[open]")
            await s.click("#newproj"); await pg.keyboard.type("Sickle-cell coursework", delay=70)
            await s.click("#savego")
            await pg.wait_for_selector(".toast", timeout=10000)
            await s.wait(1600)
            await s.scroll_to("#structcard", "start", 1200)
            await s.cap("AlphaFold 3D structure you can rotate", "coloured by AlphaFold's confidence")
            await s.drag("#viewer", 220, 40)
            await s.wait(800)
            await s.click("#color"); await s.wait(1200)
            await s.cap("Helices and strands read from the 3D model", "plus disordered regions")
            await s.drag("#viewer", -180, 60)
            await s.wait(1200)
            await s.shot("02_structure")
            await pg.wait_for_selector("svg.network", timeout=60000)
            await s.cap("Who it works with: interaction network (STRING)", "and the biological pathways it belongs to (Reactome)")
            await s.scroll_to("#ppi", "start", 1300)
            await s.move_to("svg.network .node[data-i='2']"); await s.wait(1800)
            await s.shot("03_network")

        # 2b. Projects: saved work, backup, report
        async with Scene(browser, "02b_projects") as s:
            pg = s.page
            await pg.goto(BASE + "/projects")
            await pg.wait_for_selector("#list [data-id]", timeout=30000)
            await s.click("#list [data-id]")
            await pg.wait_for_selector("#detail .item", timeout=30000)
            await s.cap("Your work is saved", "export a report or PDF with methods, back everything up in one file")
            await s.move_to("#detail a[href$='/report']"); await s.wait(1500)
            await s.shot("03b_projects")
            await s.move_to("a[href='/api/projects/backup']"); await s.wait(1800)

        # 3. Mutations
        async with Scene(browser, "03_mutations") as s:
            pg = s.page
            await pg.goto(BASE + "/mutate")
            await s.cap("Explore disease mutations", "example: the single change that causes sickle-cell anaemia")
            await s.click("#example")
            await pg.wait_for_selector("#result .stat", timeout=60000)
            await s.wait(400)
            await s.scroll_to("#result", "start", 1400)
            await s.cap("E7V: a charged residue becomes a sticky, water-repelling one",
                        "known in UniProt as Hb S, the sickle-cell variant")
            await s.wait(2600)
            await s.shot("04_mutation")
            await s.scroll_to("#heat", "center", 1300)
            await s.cap("AlphaMissense: predicted effect of every possible mutation", "blue = likely harmless · red = likely disease-causing")
            box = await pg.locator("#heat").bounding_box()
            await pg.mouse.move(box["x"] + 30, box["y"] + 60, steps=10)
            await pg.mouse.move(box["x"] + box["width"] * 0.7, box["y"] + 150, steps=70)
            await s.wait(1200)
            await s.shot("05_heatmap")

        # 4. Family tree
        async with Scene(browser, "04_tree") as s:
            pg = s.page
            await pg.goto(BASE + "/align?accs=P69905%20P68871%20P02042%20P02100%20P69891%20P02144")
            await pg.wait_for_selector("svg.tree", timeout=120000)
            await s.cap("Align many proteins and build a family tree", "human globins: β and δ pair up, myoglobin is the distant cousin")
            await s.wait(700)
            await pg.evaluate("window.scrollTo({top: document.querySelector('svg.tree').getBoundingClientRect().top + scrollY - 140, behavior: 'smooth'})")
            await s.wait(2600)
            await s.shot("06_tree")
            await s.scroll_to("#msa", "center", 1300)
            await s.cap("Multiple sequence alignment, coloured by chemistry")
            await s.wait(2000)

        # 5. 3D overlay
        async with Scene(browser, "05_overlay") as s:
            pg = s.page
            await pg.goto(BASE + "/structures?a=P69905&b=P02144")
            await pg.wait_for_selector("#viewer canvas", timeout=120000)
            await s.wait(800)
            await s.scroll_to("#out .stats", "start", 1000)
            await s.cap("Same shape, different sequence", "haemoglobin vs myoglobin: only 25% identical, yet TM-score 0.91")
            await s.wait(2200)
            await s.scroll_to("#viewer", "center", 1100)
            await s.drag("#viewer", 240, 30)
            await s.wait(1400)
            await s.shot("07_overlay")

        # 6. Classify + why (the slow explanation step is cut from the clip)
        async with Scene(browser, "06_classify") as s:
            pg = s.page
            await pg.goto(BASE + "/classify")
            await pg.wait_for_selector("#model option", state="attached", timeout=30000)
            await s.cap("Predict what a protein does with machine learning", "ESM-2 protein language model · 80% accurate on proteins it has never seen")
            await s.click("#t"); await pg.keyboard.type("P00533", delay=110)
            await s.click("#go")
            await pg.wait_for_selector("#out .card .stat", timeout=120000)
            await s.scroll_to("#out .card", "start", 1200)
            await s.wait(1600)
            await s.cap("…and see WHY the model decided", "which parts of the sequence it relied on")
            await s.click(".why")
            cut_start = s.now() + 1.2
            await pg.wait_for_selector(".whybox h3", timeout=300000)
            s.cuts.append([cut_start, s.now() - 0.3])
            await s.scroll_to(".whybox", "center", 1200)
            await s.wait(2600)
            await s.shot("08_explain")

        # 7. Lab tools
        async with Scene(browser, "07_lab") as s:
            pg = s.page
            await pg.goto(BASE + "/lab#translate")
            await s.cap("Lab tools", "DNA → protein, protease digests for mass spec, protein concentration")
            await s.click("#tex"); await s.wait(400)
            await s.click("#tgo")
            await pg.wait_for_selector("#tout tbody tr", timeout=60000)
            await s.scroll_to("#tout", "start", 1300)
            await s.cap("Haemoglobin β gene → its exact 147-amino-acid protein")
            await s.wait(2000)
            await s.shot("09_lab")
            await s.click("#tout tbody tr [data-to='digest']")
            await s.wait(500)
            await s.cap("Trypsin digest: peptide masses as a mass spectrometer sees them")
            await s.click("#dgo")
            await pg.wait_for_selector("#dout .stat", timeout=60000)
            await s.scroll_to("#dout", "start", 1200)
            await s.wait(2200)

        # 8. Learn
        async with Scene(browser, "08_learn") as s:
            pg = s.page
            await pg.goto(BASE + "/learn")
            await s.cap("Guided lessons with live data and quizzes", "for students and teachers")
            await s.wait(1800)
            await s.click("a[href='#mutation'].btn")
            await pg.wait_for_selector(".demo .stat", timeout=60000)
            await s.scroll_to(".quiz", "center", 1300)
            await s.click(".quiz label:nth-of-type(1) input")
            await s.wait(2600)
            await s.shot("10_lesson")

        await card(browser, "09_end", """<div style="font-size:52px;font-weight:800">Protein Analyzer</div>
          <div style="font-size:24px;margin-top:16px;line-height:1.5">Search · 3D structures · mutations · family trees · interactions<br>
          machine learning with explanations · lab tools · lessons<br>dark mode · save &amp; back up your work</div>
          <div style="font-size:18px;margin-top:30px;opacity:.85">Built with Claude Code (AI pair-programming)<br>
          Data: UniProt · AlphaFold DB · AlphaMissense · STRING · Reactome · InterPro · EMBL-EBI</div>""", 4.5)
        await browser.close()


def assemble(out_path: Path):
    ff = imageio_ffmpeg.get_ffmpeg_exe()
    names = sorted(p.stem for p in CLIPS.glob("*.webm"))
    parts = []
    for i, n in enumerate(names):
        cuts = json.loads((CLIPS / f"{n}.cuts.json").read_text()) if (CLIPS / f"{n}.cuts.json").exists() else []
        src = CLIPS / f"{n}.webm"
        # keep [0,a1],[b1,a2],...,[bn,end]
        keep, start = [], 0.0
        for a, b in cuts:
            keep.append((start, a)); start = b
        keep.append((start, None))
        for j, (a, b) in enumerate(keep):
            seg = CLIPS / f"seg_{i:02d}_{j}.mp4"
            cmd = [ff, "-y", "-loglevel", "error", "-ss", f"{a:.2f}"] + (["-to", f"{b:.2f}"] if b else []) + [
                "-i", str(src), "-vf", f"scale={W}:{H},fps=30,format=yuv420p", "-c:v", "libx264", "-preset", "slow",
                "-crf", "18", "-an", str(seg)]
            # first scene clip starts with a short white frame from page load; trim 0.4 s from every clip start
            if j == 0:
                cmd[cmd.index("-ss") + 1] = f"{a + 0.4:.2f}"
            subprocess.run(cmd, check=True)
            parts.append(seg)
    listfile = CLIPS / "concat.txt"
    listfile.write_text("".join(f"file '{p.as_posix()}'\n" for p in parts))
    subprocess.run([ff, "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", str(listfile), "-c:v", "libx264",
                    "-preset", "slow", "-crf", "20", "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(out_path)], check=True)


if __name__ == "__main__":
    shutil.rmtree(CLIPS, ignore_errors=True)
    asyncio.run(record())
    out = OUT / "protein_analyzer_showcase.mp4"
    assemble(out)
    print("wrote", out)
