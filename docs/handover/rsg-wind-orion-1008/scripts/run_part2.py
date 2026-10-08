"""Part 2 + 3 driver: writes rows2.pkl (gas, per phase) and rows3.pkl (dust expansion)."""
import pickle
import twophase as TP
tabs = {m: TP.part1_tables(m) for m in ('clamp', 'extrap')}
rows = TP.part2(tabs)
pickle.dump(rows, open(TP.OUT + 'rows2.pkl', 'wb'))
rows3 = TP.dust_rows(rows)
pickle.dump(rows3, open(TP.OUT + 'rows3.pkl', 'wb'))
print(len(rows), len(rows3))
