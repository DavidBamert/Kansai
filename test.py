import numpy as np
um = np.array([[0,0,0],
              [100,90,80],
              [100,120,110],
              [100,80,70],
              [0,0,0]])
uincrm = np.zeros(um.shape)

i = 1
for counter in um[0, 1:]:  # calculate uincrements
    uincrm[:, i] = um[:, i - 1] - um[:, i]
    if any(uincrm[:, i] < 0):
        for idx, n in enumerate(uincrm[:, i]):
            if uincrm[idx, i] < 0:
                uincrm[idx, i] = 0  # just at the depth where uincrm<0 is corrected

    i += 1

print(uincrm)