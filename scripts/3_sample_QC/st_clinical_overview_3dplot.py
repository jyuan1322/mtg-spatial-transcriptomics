from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# read clinical meta data from file
df = pd.read_csv("st_clinical_meta_data_V4_96samples_10022024.csv")
# df["nctx_temporal_factor_disp"] = df["nctx_temporal_factor"].replace(4, 3)
df["nctx_temporal_factor_disp"] = df["nctx_temporal_factor"]
# NOTE 1 nan value in nctx_temporal is not displayed
df = df[df["nctx_temporal_factor_disp"].notna()]
df["nctx_temporal_factor_disp"] = df["nctx_temporal_factor_disp"].astype(int)

fig = plt.figure()
ax = fig.add_subplot(111, projection = "3d")

ax.set_xlabel("nctx_temporal")
ax.set_ylabel("brain_stem_sn")
ax.set_zlabel("count")
# ax.set_xlim3d(0,10)
# ax.set_ylim3d(0,10)

# xpos = [2,5,8,2,5,8,2,5,8]
# ypos = [1,1,1,5,5,5,9,9,9]
# zpos = np.zeros(9)


xpos = df["nctx_temporal_factor_disp"].unique()
ypos = df["brain_stem_sn_factor"].unique()
# xpos_grid = np.repeat(xpos, len(ypos))
# array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4])
xpos_grid = np.array(xpos.tolist() * len(ypos))
# array([0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4])

# ypos_grid = np.array(ypos.tolist() * len(xpos))
# array([0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4])
ypos_grid = np.repeat(ypos, len(xpos))
# array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4])

zpos_grid = np.zeros(len(xpos_grid))

dx = np.ones(len(xpos_grid)) * 0.4
dy = np.ones(len(ypos_grid)) * 0.4

xpos_cent = xpos_grid - dx/2
ypos_cent = ypos_grid - dy/2

colors = {
    "Control": "#1a9641",
    "Prodromal": "#fdae61",
    "Case": "#d7191c"
}

for diag in ["Control", "Prodromal", "Case"]:
    df_sub = df[df["diagnosis"] == diag]
    counts = pd.DataFrame(df_sub[["nctx_temporal_factor_disp", "brain_stem_sn_factor"]].value_counts().reset_index())
    zpos_grid_sub = np.zeros(len(xpos_grid))
    for row in counts.itertuples():
        idx = row.brain_stem_sn_factor * len(ypos) + row.nctx_temporal_factor_disp
        zpos_grid_sub[idx] += row.count
    # https://stackoverflow.com/a/70196730
    # don't print top of box if area (count) is 0
    mask_dz = zpos_grid_sub == 0
    ax.bar3d(xpos_cent, ypos_cent, zpos_grid, dx, dy, zpos_grid_sub, color=colors[diag])
    zpos_grid += zpos_grid_sub
    

# dx = np.ones(9)
# dy = np.ones(9)
# dz = [np.random.random(9) for i in range(4)]  # the heights of the 4 bar sets

# _zpos = zpos   # the starting zpos for each bar
# colors = ['r', 'b', 'g', 'y']
# for i in range(4):
#     ax.bar3d(xpos, ypos, _zpos, dx, dy, dz[i], color=colors[i])
#     _zpos += dz[i]    # add the height of each bar to know where to start the next

# plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.show()
fig.savefig("st_clinical_overview_3dplot_10022024.pdf", bbox_inches='tight')
