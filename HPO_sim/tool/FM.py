import numpy as np

class mod(object):

    def __init__(self, INFO_PATH):
        INFO_FILE = open(INFO_PATH, "rb")
        
        buffer = np.fromfile(INFO_FILE, dtype = "uint32")

        offset = 0

        # 获取分组数量
        mod_size = buffer[offset]
        offset = offset + 1

        # 获取分组列表
        mods = buffer[offset : mod_size + offset]
        offset = offset + mod_size

        # 获取各分组数量
        sub_size = buffer[offset : mod_size + offset]
        offset = offset + mod_size

        # 生成分组字典
        self.mods = dict()
        for i in range(0, mod_size):
            self.mods[mods[i]] = set(buffer[offset : sub_size[i] + offset])
            offset = offset + sub_size[i]

        # 生成归属字典
        self.belones = dict()
        for i in self.mods.keys():
            for j in self.mods[i]:
                try:
                    self.belones[j].add(i)
                except:
                    self.belones[j]= set([i])