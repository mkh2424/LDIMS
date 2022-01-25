from PyQt5.QtWidgets import *
from PyQt5 import uic
from PyQt5 import QtCore
from PyQt5.QtGui import QGuiApplication

import os, sys, glob, shutil
import numpy as np
import datetime
import pandas as pd
import Modules.NFS_DNA as NFS_DNA


form_class = uic.loadUiType("GUI/MainSuiteForm.ui")[0]


class MainForm(QMainWindow, form_class):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        # 기본 검색 조건 설정
        self.standard_locus = ["AMEL", "D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX",
                          "D8S1179", "D21S11", "D18S51", "D2S441", "D19S433", "TH01",
                          "FGA", "D22S1045", "D5S818", "D13S317", "D7S820", "D10S1248",
                          "D1S1656", "D12S391", "D2S1338"]
        # DB 파일 불러오기, 클래스 변수 선언
        self.df_DB = {}
        self.df_query = pd.DataFrame({})
        self.dict_match = {}

        # GUI 날짜 기본값 설정
        self.date_file_to.setDate(QtCore.QDate.currentDate())
        self.date_file_from.setDate(QtCore.QDate.currentDate().addYears(-1))
        self.date_id_to.setDate(QtCore.QDate.currentDate())
        self.date_id_from.setDate(QtCore.QDate.currentDate().addYears(-1))
        # 탭 활성화/비활성화
        self.tabWidget.setTabEnabled(2, False)
        # 테이블 위젯 헤더 설정.
        list_columns = self.standard_locus
        self.table_id_query.setColumnCount(len(list_columns))
        self.table_id_query.setHorizontalHeaderLabels(list_columns)
        self.table_id_query.horizontalHeader().setDefaultSectionSize(53)
        self.table_id_query.setRowCount(1)
        self.table_id_query.setVerticalHeaderLabels(['입력'])
        self.table_id_query.verticalHeader().setMinimumWidth(70)
        for idx_col, loci in enumerate(self.standard_locus):
            item = QTableWidgetItem("")
            self.table_id_query.setItem(0, idx_col, item)
        self.table_result_query.setColumnCount(len(list_columns))
        self.table_result_query.setHorizontalHeaderLabels(list_columns)
        self.table_result_query.horizontalHeader().setDefaultSectionSize(47)
        self.table_result_query.verticalHeader().setMinimumWidth(70)
        list_columns = self.standard_locus + ['Score']
        self.table_result_target.setColumnCount(len(list_columns))
        self.table_result_query.setHorizontalHeaderLabels(list_columns)
        self.table_result_target.horizontalHeader().setDefaultSectionSize(47)
        self.table_result_target.verticalHeader().setMinimumWidth(70)

    # @ 내부 함수
    def filter_single_crosschecked(self, df_query, threshold_allele=0):
        # 한사람 allele 수가 0, 1 또는 2개인 경우) 필터링.
        score_query = df_query[self.standard_locus]
        score_query = score_query.applymap(lambda x: len(str(x).split('-')) if str(x) != '' else 0)
        score_query = score_query.applymap(lambda x: x if x in list(range(threshold_allele,3)) else np.nan)
        idx_passed = score_query.dropna().index
        return df_query.loc[idx_passed]

    def filter_overlap(self, df_query, df_target, text_progress):
        # DB 업데이트 시 같은 샘플 이름이 있는지 확인하고 중복된 것은 최신의 것으로 업데이트한 후 DB 반환
        df_sum = pd.concat([df_target, df_query])
        duplicated = df_sum.duplicated(['Sample Name'])
        if sum(duplicated) > 0:
            list_duplicated = df_sum[duplicated]['Sample Name'].unique()
            text_progress.append(f'중복되는 Sample Name이 발견되었습니다. 기존의 데이터를 제거합니다. : \n{list_duplicated}')
            df_sum = pd.concat([df_target, df_query, df_query])
            df_target_filtered = df_sum.drop_duplicates(['Sample Name'], keep=False)
        else:
            df_target_filtered = df_target
        return df_target_filtered

    def filter_date(self, df_target, date_from, date_to):
        # 지정한 기간에 저장된 DB 자료를 필터링
        df_date_filtered = df_target.copy()
        df_date_filtered['datetime'] = df_date_filtered['Date'].apply(lambda x: datetime.datetime.strptime(x, '%Y%m%d'))
        df_date_filtered = df_date_filtered[df_date_filtered['datetime'] > date_from]
        df_date_filtered = df_date_filtered[df_date_filtered['datetime'] < date_to]
        df_date_filtered.drop(['datetime'], axis=1, inplace=True)
        return df_date_filtered

    def add_query(self, df_query, df_target, text_progress):
        # 쿼리 데이터에서 교차 검증 단일 데이터와 중복 데이터를 처리하고 기존 DB에 더해서 반환
        df_query_filterd = self.filter_single_crosschecked(df_query, threshold_allele=1)
        df_target_filterd = self.filter_overlap(df_query_filterd, df_target, text_progress)
        df_sum = pd.concat([df_query_filterd, df_target_filterd], axis=0)
        df_sum.reset_index(inplace=True, drop=True)
        return df_sum

    def score_comparison(self, df_query, df_target, threshold=18):
        # 일치 여부 계산을 위해 점수표 생성하고 일치 여부를 확인한 결과인 query-index->target->indexs 쌍을 저장한 딕셔너리 반환
        dict_match = {}
        df_target_scored = df_target.copy()
        df_target_scored['score'] = 0
        for row in df_query.itertuples():
            for marker in self.standard_locus:
                df_target_scored['score'] = df_target_scored['score'] + df_target_scored[marker].apply(
                    lambda x: 1 if x == getattr(row, marker) else 0)
            df_match = df_target_scored[df_target_scored['score'] >= threshold]
            if len(df_match) > 0:
                dict_match[getattr(row, 'Index')] = df_match.index
            df_target_scored['score'] = 0
        return dict_match

    def search_DB(self, df_query, df_target, date_from, date_to, threshold=18):
        # 데이터를 조건에 맞게 필터링 하고 일치 여부 연산의 결과인 query-index->target->indexs 쌍을 반환
        df_query_filtered = self.filter_single_crosschecked(df_query)
        if len(df_query_filtered) == 0:
            QMessageBox.information(self, 'Error', 'No appropriate single profile existed.')
            return {}
        df_target_filtered = self.filter_date(df_target, date_from, date_to)
        df_target_filtered = self.filter_overlap(df_query_filtered, df_target_filtered, text_progress = self.text_result_summary)
        dict_match = self.score_comparison(df_query_filtered, df_target_filtered, threshold)
        return dict_match

    def update_table_query(self, idx_query_selected):
        """
        입력받은 index의 쿼리 데이터 프레임의 내용을 쿼리 테이블에 반영한다.
        """

        list_columns = self.standard_locus
        self.table_result_query.setHorizontalHeaderLabels(list_columns)
        self.table_result_query.setRowCount(1)
        self.table_result_query.setVerticalHeaderLabels([self.df_query.loc[idx_query_selected]['Sample Name']])
        for col_index, column in enumerate(list_columns):
            item = QTableWidgetItem(self.df_query.loc[idx_query_selected][column])
            self.table_result_query.setItem(0, col_index, item)

    def update_table_target(self, idx_query_selected, idxs_target_selected):
        """
        입력받은 index의  타겟 데이터 프레임의 내용을  타겟 테이블에 반영한다.
        """

        list_columns = self.standard_locus + ['Score']
        self.table_result_target.setHorizontalHeaderLabels(list_columns)
        self.table_result_target.setRowCount(len(idxs_target_selected))
        self.table_result_target.setVerticalHeaderLabels(self.df_DB.loc[idxs_target_selected]['Sample Name'])
        score = 0
        for row_index, idx_target_selected in enumerate(idxs_target_selected):
            for col_index, column in enumerate(list_columns):
                if column != 'Score':
                    item = QTableWidgetItem(self.df_DB.loc[idx_target_selected][column])
                    if self.df_DB.loc[idx_target_selected][column] == self.df_query.loc[idx_query_selected][column]:
                        color = QtCore.Qt.green
                        score = score + 1
                    else:
                        color = QtCore.Qt.red
                    item.setBackground(color)
                else:  # 유사도 칼럼 입력
                    item = QTableWidgetItem(str(np.round(score / 21, 1)))
                self.table_result_target.setItem(row_index, col_index, item)
            score = 0

    # @ Pyqt 이벤트
    def click_tab_resize(self, num_tab):
        """
        클릭한 탭이 Report 탭이면 프로그램 창의 크기를 키우고 다른 탭을 누르면 창 크기를 원상복구한다.

        Parameter
        ---------
        num_tab: int
            클릭한 탭의 번호(순서)
        """

        if num_tab==0:  # 파일 검색
            self.resize(654, 398)
            #self.tabWidget.resize(1535, 881)
            #self.groupBox_8.resize(1505, 721)
        elif num_tab==1: # 개별 검색
            self.resize(1255, 398)
        elif num_tab==2: # 결과
            self.resize(1437, 758)
        elif num_tab==3:
            self.resize(654, 398)

    def click_item_list_result(self, item):
        """
        리스트에서 선택된 증거물명을 쿼리로 하여 일치 결과를 각각의 테이블에 반영한다.
        """

        if item == None: return
        self.table_result_query.clear()
        self.table_result_target.clear()
        SN_query = item if isinstance(item, str) else item.text()
        idx_query_selected = self.df_query[self.df_query['Sample Name'] == SN_query].index.values[0]
        self.update_table_query(idx_query_selected)
        self.update_table_target(idx_query_selected, self.dict_match[idx_query_selected])

    def click_btn_file_tomato(self):
        """
        파일 검색 탭의 Load 버튼 이벤트. 토마토 파일을 불러오고 파싱하여 쿼리 데이터프레임에 저장.
        """
        filename_import = QFileDialog.getOpenFileName(self, "토마토 파일을 선택하세요.", os.getcwd(), 'xlsm(*.xlsm)')[0]
        if len(filename_import) == 0:
            QMessageBox.information(self, "오류", "선택된 파일이 없습니다.")
            return
        self.line_file_location.setText(filename_import)
        try:
            self.df_query = NFS_DNA.CombinedResult()
            self.df_query.load_tomato(self.line_file_location.text())
            self.df_query = self.df_query.df_profiles.astype(str)
        except:
            QMessageBox.information(self, "Error", f'Failed to load the tomato file - Filename : {self.line_file_location.text()}')
            return 0
        QMessageBox.information(self, "Load", 'Load complete.')

    def click_btn_file_id(self):
        """
        개별 검색 탭의 검색 버튼의 클릭 이벤트. 쿼리와 타겟 데이터프레임을 생성하고 결과를 출력한다.
        """
        self.df_DB = pd.read_csv('local_DB.csv', header=0, index_col=0, dtype=str)
        self.list_result_matched.clear()
        self.table_result_query.clear()
        self.table_result_target.clear()
        self.text_result_summary.clear()
        # table_id_query로 부터 df_query 생성
        df_query = {'Sample Name': ['Query'], 'Date': [datetime.date.today().strftime('%Y%m%d')]}
        for idx_col, loci in enumerate(self.standard_locus):
            df_query[loci] = self.table_id_query.item(0,idx_col).text()
        self.df_query = pd.DataFrame.from_dict(df_query)
        # 검색일에 해당하는 DB 데이터만 추출하고 서치
        date_from = datetime.datetime.strptime(self.date_id_from.date().toString('yyyyMMdd'), '%Y%m%d')
        date_to = datetime.datetime.strptime(self.date_id_to.date().toString('yyyyMMdd'), '%Y%m%d')
        self.dict_match = self.search_DB(self.df_query, self.df_DB, date_from, date_to,
                                         threshold=self.sb_id_threshold.value())
        if len(self.dict_match.keys()) == 0:
            QMessageBox.information(self, '결과', '일치 건이 없습니다.')
            return 0
        # 결과 텍스트 창에 검색 결과 요약 출력
        self.text_result_summary.append(f'{len(self.dict_match.keys())} 건의 일치 건이 있습니다.\n-------------------------------\n')
        for idx_query, idxs_match in self.dict_match.items():
            for idx_match in idxs_match:
                sn_query = self.df_query.loc[idx_query]['Sample Name']
                sn_match = self.df_DB.loc[idx_match]['Sample Name']
                self.text_result_summary.append(f'{sn_query} => {sn_match}')
        self.text_result_summary.append('--------------------------------\n')
        # 일치한 샘플 이름을 리스트에 반영하고 좌위 테이블 출력
        if len(self.dict_match.keys())>0:
            self.list_result_matched.addItems(self.df_query.loc[self.dict_match.keys()]['Sample Name'])
            self.click_item_list_result(self.df_query.loc[self.dict_match.keys()]['Sample Name'].values[0])
        self.tabWidget.setTabEnabled(2, True)
        self.tabWidget.setCurrentIndex(2)

    def click_btn_file_search(self):
        """
        파일 검색 탭의 검색 버튼의 클릭 이벤트. 쿼리와 타겟 데이터프레임을 생성하고 결과를 출력한다.
        """
        self.df_DB = pd.read_csv('local_DB.csv', header=0, index_col=0, dtype=str)
        self.list_result_matched.clear()
        self.table_result_query.clear()
        self.table_result_target.clear()
        self.text_result_summary.clear()
        # 검색일에 해당하는 DB 데이터만 추출하고 서치
        date_from = datetime.datetime.strptime(self.date_file_from.date().toString('yyyyMMdd'), '%Y%m%d')
        date_to = datetime.datetime.strptime(self.date_file_to.date().toString('yyyyMMdd'), '%Y%m%d')
        self.dict_match = self.search_DB(self.df_query, self.df_DB, date_from, date_to, threshold=self.sb_file_threshold.value())
        self.text_result_summary.append(f'{len(self.dict_match.keys())} 건의 일치 건이 있습니다.\n-------------------------------')
        if len(self.dict_match.keys()) != 0:
            # 결과 텍스트 창에 검색 결과 요약 출력
            for idx_query, idxs_match in self.dict_match.items():
                for idx_match in idxs_match:
                    sn_query = self.df_query.loc[idx_query]['Sample Name']
                    sn_match = self.df_DB.loc[idx_match]['Sample Name']
                    self.text_result_summary.append(f'{sn_query} => {sn_match}')
            self.text_result_summary.append('--------------------------------')
            # 일치한 샘플 이름을 리스트에 반영하고 좌위 테이블 출력
            self.list_result_matched.addItems(self.df_query.loc[self.dict_match.keys()]['Sample Name'])
            self.click_item_list_result(self.df_query.loc[self.dict_match.keys()]['Sample Name'].values[0])
        else:
            QMessageBox.information(self, '결과', '일치 건이 없습니다.')
        self.tabWidget.setTabEnabled(2, True)
        self.tabWidget.setCurrentIndex(2)
        # 저장에 체크되어 있으면 쿼리 데이터를 기존 DB에 저장.
        if self.cb_file_db_save.isChecked():
            time_now = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
            filename_backup = f'backup_{time_now}.csv'
            shutil.copy('local_DB.csv', './Backup/' + filename_backup)
            self.text_result_summary.append(f'기존 DB를 백업합니다. : ./Backup/{filename_backup}')
            self.text_result_summary.append('--------------------------------')
            df_sum = self.add_query(self.df_query, self.df_DB, text_progress=self.text_result_summary)
            df_sum.sort_values(by=['Date'], axis=0, ascending=True, inplace=True)
            df_sum.reset_index(drop=True, inplace=True)
            df_sum.to_csv("local_DB.csv", sep=',')
            self.text_result_summary.append('저장 완료.')

    def update_files(self, files_update):
        time_now = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        filename_backup = f'backup_{time_now}.csv'
        shutil.copy('local_DB.csv', './Backup/' + filename_backup)
        self.text_tool_progress.append(f'기존 DB를 백업합니다. : ./Backup/{filename_backup}')
        self.text_tool_progress.append('--------------------------------------------------------------------------------------------\n')
        self.text_tool_progress.append(f'총 {len(files_update)}의 파일을 업데이트 합니다.')

        df_DB = pd.read_csv('local_DB.csv', header=0, index_col=0, dtype=str)
        cnt_db_old = len(df_DB)
        cnt_stored = 0
        for idx, filename in enumerate(files_update):
            QGuiApplication.processEvents()
            try:
                query = NFS_DNA.CombinedResult()
                query.load_tomato(filename)
            except:
                message = f'{idx + 1} : 다음 파일을 불러오는 데 실패했습니다. - 파일명 : {filename}'
                self.text_tool_progress.append(message)
            else:
                cnt_stored += 1
                df_DB = self.add_query(query.df_profiles, df_DB, text_progress=self.text_tool_progress)
                message = f'{idx + 1} : {filename}이 저장되었습니다.'
                self.text_tool_progress.append(message)
        df_DB.sort_values(by=['Date'], axis=0, ascending=True, inplace=True)
        df_DB.reset_index(drop=True, inplace=True)
        df_DB.to_csv("local_DB.csv", sep=',')
        self.text_tool_progress.append('--------------------------------------------------------------------------------------------\n')
        message = f'총 {len(files_update)}개의 파일 중 {cnt_stored}개의 파일을 업데이트했습니다.'
        self.text_tool_progress.append(message)
        message = f'{len(df_DB) - cnt_db_old}개의 데이터가 추가되었습니다.\nDB에 저장된 프로필 수는 {len(df_DB)}입니다.'
        self.text_tool_progress.append(message)

    def click_btn_tool_update_from_files(self):
        self.text_tool_progress.clear()
        files_update = QFileDialog.getOpenFileNames(self, '토마토 파일을 선택하세요.',os.getcwd(), 'xlsm(*.xlsm)')[0]
        if len(files_update)==0:
            QMessageBox.information(self, "오류", "선택된 파일이 없습니다..")
            return
        self.update_files(files_update)

    def click_btn_tool_update_from_folder(self):
        self.text_tool_progress.clear()
        folder_update = QFileDialog.getExistingDirectory(self, '파일이 저장된 폴더를 선택하세요.', os.getcwd())
        files_update = glob.glob(folder_update+'/*.xlsm')
        if len(files_update)==0:
            QMessageBox.information(self, "오류", "선택된 파일이 없습니다..")
            return
        self.update_files(files_update)

    def click_btn_tool_load_backup(self):
        self.text_tool_progress.clear()
        file_backup = QFileDialog.getOpenFileName(self, '백업 파일을 선택하세요.', os.getcwd()+'/Backup/', 'csv(*.csv)')[0]
        if len(file_backup)==0:
            QMessageBox.information(self, "오류", "선택된 파일이 없습니다..")
            return
        time_now = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        filename_backup = f'backup_{time_now}.csv'
        shutil.copy('local_DB.csv', './Backup/' + filename_backup)
        self.text_tool_progress.append(f'기존 DB를 백업합니다. : ./Backup/{filename_backup}')
        self.text_tool_progress.append('--------------------------------------------------------------------------------------------\n')
        df_DB = pd.read_csv(file_backup, header=0, index_col=0, dtype=str)
        df_DB.to_csv("local_DB.csv", sep=',')
        message = f'다음 파일로 DB 파일을 교체했습니다 : {file_backup}'
        self.text_tool_progress.append(message)

    def click_btn_id_load_clipboard(self):
        try:
            datas = pd.read_clipboard(sep='\n').columns[0].split('\t')
            range_copy = range(min([len(datas),21]))
            for idx_col in range_copy:
                item = QTableWidgetItem(datas[idx_col])
                self.table_id_query.setItem(0, idx_col, item)
        except:
            QMessageBox.information(self, '오류', '올바르지 않은 데이터입니다.')

def except_hook(cls, exception, traceback): # for PyQt5.5 debugging
    sys.__excepthook__(cls, exception, traceback)

if __name__ == "__main__":
        app = QApplication(sys.argv)
        entry_GUI = MainForm()
        #entry_GUI.setFixedSize(entry_GUI.size())  # 창 크기 고정
        entry_GUI.move(10, 10)
        entry_GUI.show()
        sys.excepthook = except_hook    # for PyQt5.5 debugging
        app.exec_()
