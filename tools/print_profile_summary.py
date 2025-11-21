import pstats, io, json, sys, os

if len(sys.argv) < 2:
    print(json.dumps({"error":"usage: print_profile_summary.py <profile_file>"}))
    sys.exit(1)

path = sys.argv[1]
if not os.path.exists(path):
    print(json.dumps({"error":"profile file not found","path":path}))
    sys.exit(2)

st = pstats.Stats(path)
st.strip_dirs()

entries = []
for (filename, line, func), (cc, nc, tt, ct, callers) in st.stats.items():
    entries.append({
        "func": func,
        "file": filename,
        "line": line,
        "cc": cc,
        "nc": nc,
        "tottime": tt,
        "cumtime": ct,
    })

by_cum = sorted(entries, key=lambda e: e["cumtime"], reverse=True)[:30]
by_tot = sorted(entries, key=lambda e: e["tottime"], reverse=True)[:30]

sio = io.StringIO(); p = pstats.Stats(path, stream=sio); p.strip_dirs(); p.sort_stats('cumulative').print_stats(25)
hr_cum = sio.getvalue()
sio = io.StringIO(); p = pstats.Stats(path, stream=sio); p.strip_dirs(); p.sort_stats('tottime').print_stats(25)
hr_tot = sio.getvalue()

print(json.dumps({"top_cumulative": by_cum, "top_tottime": by_tot, "hr_cumulative": hr_cum, "hr_tottime": hr_tot}))
