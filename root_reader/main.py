import os, sys

from collections import OrderedDict

import ROOT


class ROOTreader(object):
    def __init__(self, input_file):
        self.file = ROOT.TFile(input_file)
        self.root_files_path = '../assets/root_files/'
        self.lm_files_path = '../assets/lm_files/'

    def events(self, branch="Coincidences"):
        """
        Take all event ID from root file
        :param branch:
        :return:
        """

        t = self.file.Get(branch)
        event_list = []
        for event in t:
            if event.eventID1 == event.eventID2:
                event_list.append(event.eventID1)
            else:
                print "False con %d %d" % (event.eventID1, event.eventID2)

        return set(event_list)


    def time_steps(self, event_list=None, branch="Coincidences", window=3 * 10 ** (-9)):
        """
        Take list od first time from concidenses from event id
        :param event_list:
        :param branch:
        :param window:
        :return:
        """
        t = self.file.Get(branch)
        time_list = []
        if event_list:
            for event in t:
                if abs(event.time1 - event.time2) <= window:
                    if event.eventID1 in event_list:
                        time_list.append(min(event.time1, event.time2))
                else:
                    print abs(event.time1 - event.time2)
                    pass
                    # print "bad time window %d %d" % (event.time1, event.time2)
        else:
            for event in t:
                if abs(event.time1 - event.time2) <= window:
                    time_list.append(min(event.time1, event.time2))
                else:
                    pass
                    print "Bad time window %d" % (event.time1, event.time2)
        return time_list


    def hist(self, time_list, widow=3 * 10 ** (-9), branch="Singles"):
        """
        Draw a histogram of energy in branch in time window.

        :param time_list: list of coincidence time
        :param widow: time in s
        :param branch:
        :return:
        """
        c1 = ROOT.TCanvas('c1', 'Dynamic Filling Example', 200, 10, 700, 500)
        c1.SetFillColor(42)
        c1.GetFrame().SetFillColor(21)
        c1.GetFrame().SetBorderSize(6)
        c1.GetFrame().SetBorderMode(-1)
        t = self.file.Get(branch)
        print t.ls()
        query = ''
        for time in time_list:
            time_w = time + widow
            query = '(time>=%.14f) && (time<=%.14f)' % (time-widow, time_w)
            print  query
            if branch == 'Hits':
                t.Draw('edep', query)
            else:
                t.Draw('energy', query)
            c1.Modified()
            c1.Update()
            raw_input("Press enter to continue...")


    def root2lm(self, branch="Coincidences", out_put="../assets/lm_files/out_put.txt", event_list=[]):
        """
        Parse Coincidences to lm

        :param branch:
        :param out_put:
        :param event_list: parse only slected events
        :return:
        """
        list = []
        t = self.file.Get(branch)
        if not event_list:
            for event in t:
                list.append("%.2f %.2f %.2f %.5f %.2f %.2f %.2f %.5f %d %d %.5f %.5f %d %d " % (
                    event.globalPosX1 / 10., event.globalPosY1 / 10., event.globalPosZ1 / 10., event.time1 * 1e12,
                    event.globalPosX2 / 10., event.globalPosY2 / 10., event.globalPosZ2 / 10., event.time2 * 1e12,
                    event.rsectorID1, event.rsectorID2, event.energy1 * 1000, event.energy2 * 1000, event.eventID1,
                    event.eventID2))
        else:
            for event in t:
                if event.eventID1 in event_list:
                    list.append("%.2f %.2f %.2f %.5f %.2f %.2f %.2f %.5f %d %d %.5f %.5f %d %d " % (
                        event.globalPosX1 / 10., event.globalPosY1 / 10., event.globalPosZ1 / 10., event.time1 * 1e12,
                        event.globalPosX2 / 10., event.globalPosY2 / 10., event.globalPosZ2 / 10., event.time2 * 1e12,
                        event.rsectorID1, event.rsectorID2, event.energy1 * 1000, event.energy2 * 1000, event.eventID1,
                        event.eventID2))
        out_file = open(out_put, 'w')
        for line in list:
            out_file.write("%s\n" % line)


    def Hits2ODict(self, trashhold=0.01):
        """
        Return hits in OrderedDict
        :param trashhold:
        :return:
        """


        Hits = self.file.Get('Hits')
        ordered_hits = OrderedDict()
        for Hit in Hits:
            # print Hit.processName,Hit.processName, Hit.nPhantomRayleigh, Hit.nCrystalRayleigh, Hit.PDGEncoding
            #
            # if Hit.eventID > 100:
            #     break
            if Hit.edep >= trashhold and "compt" in Hit.processName and Hit.nPhantomRayleigh == 0 and\
                            Hit.nCrystalRayleigh == 0 and Hit.PDGEncoding == 22:
                if ordered_hits.has_key(Hit.eventID):
                    ordered_hits[Hit.eventID].append({"edep": Hit.edep, "time": Hit.time})
                else:
                    ordered_hits[Hit.eventID] = [{"edep": Hit.edep, "time": Hit.time}]
        return ordered_hits

    def EventpruneCoincidences_ODict(self,ordered_hits, trashhold=0.2, con_branch="Coincidences", limit=3,window=3 * 10 ** -9):
        """
        Return list of event J-PET Coincidences

        :param ordered_hits:
        :param trashhold:
        :param con_branch:
        :param limit:
        :param window:
        :return:
        """

        Coincidences = self.file.Get(con_branch)
        event_list = []
        bad_list = []
        for Coincidence in Coincidences:
            event = Coincidence.eventID1
            if ordered_hits.get(event):
                time = min(map(lambda x: x['time'], ordered_hits.get(event)))
            p = 0
            trash = 0
            for ev in range(event - limit, event + limit):
                event_ = ordered_hits.get(ev)
                if event_:
                    sorted(event_, key=lambda x:x['time'])
                    for hit in event_:
                        if abs(hit['time']-time) <= window:
                            time = hit['time']
                            p += 1
                            if hit['edep'] >= trashhold:
                                trash += 1

            if p > 3 or trash != 2:
                bad_list.append(event)
            else:
                event_list.append(event)
        print 'pruned coincidences', len(bad_list)
        return event_list

    def show_element(self,event, branch ='Hits'):
        elements = self.file.Get(branch)
        for element in elements:
            if event == element.eventID:
                print element.Show()

    def prunedTree(self,events, output, branch_list):
        ofile = ROOT.TFile(os.path.join(self.root_files_path, output), 'recreate')
        for branch in branch_list:
            tree = self.file.Get(branch)
            print 'events len: ', len(events)


            new_tree = tree.CloneTree(0)
            for event in tree:
                if "Coincidences" in branch:
                    if event.eventID1 in list(events):
                        new_tree.Fill()
                else:
                    if event.eventID in list(events):
                        new_tree.Fill()
            print new_tree.Print()

            ofile.cd()
            new_tree.Write()
        print ofile.ls()
        ofile.Close()



class LMreader(object):
    def __init__(self, input_file):
        path_file = os.path.dirname(os.path.realpath(__file__))
        try:
            self.file = open(os.path.join(path_file, input_file), 'r')
        except IOError:
            print("No such file as %s" % os.path.join(path_file, input_file))
            sys.exit(2)


    def events(self):
        """
        Return set of events
        :return:
        """
        event_list = []
        for line in self.file.readlines():
            line = list(map(float, line.split()))
            try:
                if line[13] == line[12]:
                    event_list.append(int(line[12]))
                else:
                    print "False con %d %d" % (line[13], line[12])

            except IndexError:
                print("Bad structure. X1 Y1 Z1 time1 X2 Y2 Z2 time2 ID1 ID2 Energy1 Energy2 Type")
                sys.exit(2)
        return set(event_list)


# def show(lis):
#     time = min(map(lambda x: x['time'], lis))
#     return time, map(lambda x: {'edep': '%.3f' % x['edep'], 'time': x['time'] - time}, lis)
#
# def order(time, number, order_dict):
#     i = number + 1
#     j = number - 1
#     while order_dict.get(i) == None:
#         i += 1
#     while order_dict.get(j) == None:
#         j -= 1
#     print i, map(lambda x: {'edep': '%.3f' % x['edep'], 'time': x['time'] - time}, order_dict.get(i))
#     print j, map(lambda x: {'edep': '%.3f' % x['edep'], 'time': x['time'] - time}, order_dict.get(j))
#
# for e in [1516772, 1044906, 1577904, 1985237, 1762140, 1944469]:
#     time, ev = show(ordered_hits[e])
#     print e, ev
#     order(time, e, ordered_hits)

if __name__ == "__main__":

    ROOT_FILES = "../assets/root_files/"
    LM_FILES = "../assets/lm_files/"

    LM = os.path.join(LM_FILES,'con_z500.txt')

    for file_name in ['000','0020','2500','4000']:
        TREE = os.path.join(ROOT_FILES,'sym%s.root' % file_name)
        # TREE = os.path.join(ROOT_FILES,'source25.root')

        cylindrical = ROOTreader(TREE)
        # lm = LMreader(LM)
        # cylindrical.root2lm(out_put=os.path.join(LM_FILES,'25sourceHE.txt'),branch="HECoincidences")
        # e_root = cylindrical.events(branch="HECoincidences")
        # e_lm = lm.events()
        #
        # cylindrical.prunedTree(e_lm)
        #
        # print len(e_root), len(e_root.difference(e_lm))
        # print len(e_lm), len(e_lm.difference(e_root))

        # ordered_hits = cylindrical.Hits2ODict()
        # ev_coincidences = cylindrical.EventpruneCoincidences_ODict(ordered_hits,
        #                                                            con_branch="HECoincidences", limit=1)
        # cylindrical.prunedTree(ev_coincidences,'cut_sym%s.root' % file_name,['Hits',"HECoincidences" ])
        for b in ["HECoincidences"]:
            cylindrical.root2lm(out_put=os.path.join(LM_FILES, 'lm_%s%s.txt' % (file_name,b[:2])), branch=b)
    # print len(ev_coincidences), len(set(ev_coincidences).difference(e_lm))
    # print set(ev_coincidences).difference(e_lm)

