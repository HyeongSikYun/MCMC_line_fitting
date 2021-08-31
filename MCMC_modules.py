import numpy as np

def rand_norm(mean,std,num):
  pick = np.random.normal(0.0,1.0,num*10)
  choice = np.random.choice(pick,num)
  rand_numbers = (choice * std) + mean
  return rand_numbers

def rand_uniform(num):
  pick = np.random.uniform(low=0.0,high=1.0,size=num*10)
  choice = np.random.choice(pick,num)
  return choice

def judge_acceptance(p):
  test_p = rand_uniform(1)
  out = False
  if test_p <= p:
    out = True
  return out

def make_triangle_image(val1, val2):
  min_max_val1 = [np.min(val1),np.max(val1)]
  min_max_val2 = [np.min(val2),np.max(val2)]
  im_mat = np.zeros([50,50])
  val1_arr = np.linspace(min_max_val1[0],min_max_val1[1],51)
  val2_arr = np.linspace(min_max_val2[0],min_max_val2[1],51)
  for i in range(50):
    for j in range(50):
      idx1 = val1 > val1_arr[i]
      idx2 = val1 < val1_arr[i+1]
      idx3 = val2 > val2_arr[j]
      idx4 = val2 < val2_arr[j+1]
      idx12 = np.logical_and(idx1,idx2)
      idx34 = np.logical_and(idx3,idx4)
      idx = np.logical_and(idx12,idx34)
      im_mat[j,i] = np.sum(idx)
  return val1_arr, val2_arr, im_mat
