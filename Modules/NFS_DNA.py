import pandas as pd
import os.path, datetime # 파일의 수정일을 얻기 위함
import re

class STRProfile():
    """
    법의학 실험에서 생성된 STR 데이터를 저장하고 가공하는 클래스

    Attributes
    ----------
    id : str
        프로파일의 ID
    profile : dict
        좌위-좌위값을 키-값으로 가지는 딕셔너리
    __TA_THRESHOLD : int
        혼합 프로파일 판정시 Triallelic을 몇 개까지 용인할 것인지를 정하는 내부변수

    Methods
    --------
    __find_common_locus(target_profile, query_profile)
        두 STRProfile에 동일하게 포함된 좌위명을 반환하는 내부함수
    compare(query)
        입력받은 프로파일과의 동일 여부를 반환
    check_inclusion(query)
        입력받은 프로파일의 포함 여부를 반환
    check_MX()
        해당 프로파일이 혼합형인지 여부를 반환
    union_profiles(query)
        입력받은 프로파일과 해당 프로파일의 혼합형 프로파일을 생성하여 반환
    __check_special_case(self, loci, allele)
        해당 좌위의 좌위값이 데이터베이스 운영지침 상 특수 케이스인지 여부를 반환하는 내부 함수
    transform_to_str(flag_homo_duplication=True)
        profile 객체를 감정서에 들어갈 포멧으로 변환하여 반환.
        감정서 상 좌위테이블의 기타 란에 들어갈 문구가 있다면 함께 반환
    """

    def __init__(self, id="", profile={}):
        self.id = id
        self.profile = profile
        self.__TA_THRESHOLD = 0

    def __find_common_locus(self, target_profile, query_profile):
        """
        두 STRProfile에 동일하게 포함된 좌위명을 반환
        """
        locus_target = set(target_profile.keys())
        locus_query = set(query_profile.keys())
        return locus_target.intersection(locus_query)

    def rename(self, new_id):
        old_id = self.id
        self.id = new_id
        print(f"ID CHANGED : {old_id} -> {new_id}")

    def input_loci(self, loci, alleles):
        self.profile[loci] = alleles

    def input_locus(self, new_profile):
        self.profile = {**self.profile, **new_profile}

    def delete_loci(self, loci):
        del self.profile[loci]

    def delete_locus(self, locus):
        for loci in locus:
            del self.profile[loci]

    def compare(self, query):
        for loci in self.__find_common_locus(self.profile, query.profile):
            if self.profile[loci] != query.profile[loci]:
                return False
        return True

    def check_inclusion(self, query):
        for loci in self.__find_common_locus(self.profile, query.profile):
            alleles_target = set(self.profile[loci])
            alleles_query = set(query.profile[loci])
            if alleles_query.issubset(alleles_target) == False:
                return False
        return True

    def check_MX(self):
        cnt_ta = 0
        for allele in self.profile.values():
            if len(allele) > 2:
                cnt_ta = cnt_ta + 1
        if cnt_ta > self.__TA_THRESHOLD:
            return True
        else:
            return False

    def union_profiles(self, query):
        union_profile = STRProfile()
        union_profile.rename(f"{self.id} + {query.id}")
        for loci in self.__find_common_locus(self.profile, query.profile):
            alleles_target = set(self.profile[loci])
            alleles_query = set(query.profile[loci])
            union_profile.profile[loci] = list(alleles_target.union(alleles_query))
        return union_profile

    def __check_special_case(self, loci, allele):
        decimal_place = allele.split('.')[1]
        if decimal_place == "2":
            return True
        if loci == "TH01" and allele == "9.3":
            return True
        elif loci == "D2S441" and allele == "9.1":
            return True
        elif loci == "D1S1656" and (allele == "17.3" or allele == "18.3"):
            return True
        elif loci == "Penta E" and (decimal_place == "2" or decimal_place == "3"):
            return True
        elif loci == "Penta D" and (decimal_place == "2" or decimal_place == "3"):
            return True
        else:
            return False

    def transform_to_str(self, flag_homo_duplication=True):
        string_profile = {}
        flag_NC = False
        flag_ND = False
        cnt_microvariant = 0
        str_etc = []
        str_etc_microvariant = []
        temp_profile={}
        for loci, alleles in self.profile.items():
            temp_alleles = []
            if alleles[0]=='ND':
                temp_alleles.append('ND')
                flag_ND=True
            elif alleles[0]=='NC':
                temp_alleles.append('NC')
                flag_NC=True
            else:
                for allele in alleles:
                    if allele.find("OL")!=-1:
                        continue
                    modified_allele = allele
                    if allele.find('.')!=-1:
                        if not self.__check_special_case(loci, allele):#특별 케이스가 아니면
                            cnt_microvariant = cnt_microvariant + 1
                            decimal_place = allele.split('.')[1]
                            if decimal_place == "1":
                                modified_allele = str(int(float(allele)))
                            else:
                                modified_allele = str(int(float(allele)+1))
                            modified_allele = modified_allele+"*"*cnt_microvariant
                            str_etc_microvariant.append("*"*cnt_microvariant + " : 미세변이 (검출값 {0})".format(allele))
                    temp_alleles.append(modified_allele)
            temp_profile[loci]=temp_alleles
        if flag_ND==True: str_etc.append('ND: 디엔에이형이 검출되지 않음.')
        if flag_NC==True: str_etc.append('NC: 디엔에이형을 결정할 수 없음.')

        if self.check_MX()==True:
            str_etc.append('/ : 혼합 디엔에이형.')
            for loci, alleles in temp_profile.items():
                string_profile[loci]='/'.join([str(element) for element in alleles])
        else:
            for loci, alleles in temp_profile.items():
                if flag_homo_duplication==True and len(alleles)==1:
                    alleles = alleles*2
                string_profile[loci]='-'.join([str(element) for element in alleles])
        str_etc=str_etc+str_etc_microvariant
        return string_profile, '\r\n'.join(str_etc)


class CombinedResult():
    """
    Tomato 엑셀파일의 Combined_result 혹은 Genemapper Result 파일의 데이터를 STRProfile 클래스에 파싱하여
    사건번호-STRProfile을 키-값으로 가지는 딕셔너리 형태로 저장하는 클래스

    Attributes
    ----------
    profiles : dict
        사건번호-STRProfile을 키-값을로 가지는 딕셔너리
    info : DataFrame
        Tamato의 Combined_result가 가지는 좌위 값 외의 데이터를 저장하는 데이터프레임
    kit : str
        읽어올 결과에 사용된 kit 종류
    list_marker_ordered : list
        표준감정서의 좌위테이블에 좌위가 들어가는 순서를 저장한 list

    Methods
    --------
    load_tomato(self, filename)
        Tomato 엑셀 파일의 Cominbed_result 데이터를 분석을 위한 형태로 가공하여 샘플명-STRProfile 객체를
        키-밸류 값으로 가지는 딕셔너리로 만들어 저장한다.
    load_genemapper(self, filename):
        GeneMapper 결과 파일을 분석을 위한 형태로 가공하여 샘플명-STRProfile 객체를
        키-밸류 값으로 가지는 딕셔너리로 만들어 저장한다.
    """

    dict_markers = {"GF/PPF": ["AMEL", "D3S1358", "vWA", "D16S539", "CSF1PO", "TPOX",
                               "D8S1179", "D21S11", "D18S51", "D2S441", "D19S433", "TH01",
                               "FGA", "D22S1045", "D5S818", "D13S317", "D7S820", "D10S1248",
                               "D1S1656", "D12S391", "D2S1338", "Penta E", "Penta D", "SE33"],
                    "Y23": ['DYS576', 'DYS389 I', 'DYS448', 'DYS389 II', 'DYS19', 'DYS391',
                            'DYS481', 'DYS533', 'DYS438', 'DYS437', 'DYS570',
                            'DYS635', 'DYS390', 'DYS439', 'DYS392', 'DYS393',
                            'DYS458', 'DYS385', 'DYS456', 'Y GATA H4']}  # Globalfiler, Y23 기준 marker 순서, 출력폼 index 생성시 참조

    def __init__(self, kit="GF/PPF"):
        self.profiles = {}
        self.df_profiles = pd.DataFrame()
        self.info = pd.DataFrame()
        self.kit = kit
        self.list_marker_ordered = self.dict_markers[kit]

    def load_tomato(self, filename):
        """
            Tomato 엑셀 파일의 Cominbed_result 데이터를 분석을 위한 형태로 가공하여 샘플명-STRProfile 객체를
            키-밸류 값으로 가지는 딕셔너리로 만들어 저장한다.

            Parameters
            ----------
            filename : str
                읽어들일 GeneMapper 결과 파일
        """
        if self.kit=="GF/PPF":
            df_tomato = pd.read_excel(filename, sheet_name="CombinedResult", header=1)
            # cross-check된 결과(Sample Id가 공란)만 필터
            df_crosschecked = df_tomato[df_tomato['Sample ID'].isna()]
            # Sample Name이 사건번호의 포멧에 일치하는 데이터만 남김
            p = re.compile('\d+[-]\w[-]\d+')
            cond1 = df_crosschecked['Sample Name'].apply(lambda x: True if p.match(x) else False)
            df_crosschecked = df_crosschecked[cond1]
            # 중복되는 Sample Name 제거
            df_crosschecked.drop_duplicates(['Sample Name'], keep='first', inplace=True)
            # 칼럼명 Amelogenin->AMEL 변경
            df_crosschecked.rename({'Amelogenin': 'AMEL'}, axis='columns', inplace=True)
            # 좌위 추출 및 편집
            df_locus = df_crosschecked.loc[:, ['Sample Name'] + self.list_marker_ordered]
            df_locus.fillna("", inplace=True)
            # 데이터프레임 형태로 데이터를 저장
            self.df_profiles = df_locus.copy()
            mtime = datetime.datetime.fromtimestamp(os.path.getmtime(filename))
            self.df_profiles['Date'] = mtime.strftime('%Y%m%d')
            # STR profile 객체의 딕셔너리로 데이터를 저장
            df_locus.set_index('Sample Name', inplace=True)
            df_locus = df_locus.applymap(lambda x: str(x).split('-'))
            dict_temp = df_locus.to_dict(orient='index')
            for sample_name in dict_temp.keys():
                self.profiles[sample_name] = STRProfile(id=sample_name, profile=dict_temp[sample_name])
            # 기타 정보 저장
            df_info = df_crosschecked.loc[:, ['Sample Name', 'DB Type 1', 'DB Type 2', 'Matching Probability']]
            df_info.set_index('Sample Name', inplace=True)
            self.info = df_info
            print(f"{filename} loaded.")
        elif self.kit=="Y23":
            df_tomato = pd.read_excel(filename, sheet_name="CombinedResult", header=1)
            # Sample Name이 사건번호의 포멧에 일치하는 데이터만 남김
            p = re.compile('\d+[-]\w[-]\d+')
            cond1 = df_tomato['Sample Name'].apply(lambda x: True if p.match(x) else False)
            df_tomato = df_tomato[cond1]
            # 좌위 추출 및 편집
            print(set(df_tomato.columns))
            print(set(self.list_marker_ordered))

            df_locus = df_tomato.loc[:, ['Sample Name'] + self.list_marker_ordered]
            df_locus.fillna("", inplace=True)
            df_locus.set_index('Sample Name', inplace=True)
            df_locus = df_locus.applymap(lambda x: str(x).split('-'))
            dict_temp = df_locus.to_dict(orient='index')
            for sample_name in dict_temp.keys():
                self.profiles[sample_name] = STRProfile(id=sample_name, profile=dict_temp[sample_name])
            print(f"{filename} loaded.")

    def load_genemapper(self, filename):
        """
            GeneMapper 결과 파일을 분석을 위한 형태로 가공하여 샘플명-STRProfile 객체를
            키-밸류 값으로 가지는 딕셔너리로 만들어 저장한다.

            Parameters
            ----------
            filename : str
                읽어들일 GeneMapper 결과 파일
        """
        df = pd.read_csv(filename, sep='\t', dtype=str, engine='python')
        df = df.filter(regex=r'Sample Name|Marker|Allele', axis=1)  # 필요한 column만 추출
        df = df[df['Marker'].isin(self.list_marker_ordered)]  # 필요한 Marker만 추출
        df['ProcessedAllele'] = df.filter(regex=r'Allele', axis=1).apply(lambda x: x.dropna().values.tolist(),
                                                                         axis=1)  # allele1, allele2... 식으로 되어있는 allele 값을 모아서 하나의 list로 만들어 저장
        df = df[['Sample Name', 'Marker', 'ProcessedAllele']]
        # Sample Name이 사건번호의 포멧에 일치하는 데이터만 남김
        p = re.compile('\d+[-]\w[-]\d+')
        cond1 = df['Sample Name'].apply(lambda x: True if p.match(x) else False)
        df = df[cond1]
        # Marker를 칼럼으로, Sample Name을 index로 pivot
        df = df.pivot(index='Sample Name', columns='Marker', values='ProcessedAllele')
        # 데이터프레임 형태로 데이터를 저장

        self.df_profiles = df.copy()
        self.df_profiles = self.df_profiles.applymap(lambda x: '-'.join(x))
        self.df_profiles = self.df_profiles.reset_index()
        mtime = datetime.datetime.fromtimestamp(os.path.getmtime(filename))
        self.df_profiles['Date'] = mtime.strftime('%Y%m%d')
        # STR profile 객체의 딕셔너리로 데이터를 저장
        dict_temp = df.to_dict(orient='index')
        for sample_name in dict_temp.keys():
            self.profiles[sample_name] = STRProfile(id=sample_name, profile=dict_temp[sample_name])
        print(f"{filename} loaded.")