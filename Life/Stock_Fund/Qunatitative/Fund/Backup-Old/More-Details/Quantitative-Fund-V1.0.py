import sys
reload(sys)
sys.setdefaultencoding('utf-8')
                       
def html_download(self, url):
    if url is None:
        return None
    
    try:
        response = requests.get(url)
    except Exception as e:
        print("Open the failed, error :{}".format(e))
        return None
    
    if response.status_code != 200:
        return None
    
    return response.content

def html_extract_content(self,html_cont):
    funds_text=[]
    if html_cont is None:
        return
    soup=BeautifulSoup(html_cont,'html.parser',from_encoding='gb18030')
    #Get all the fund ID.
    print(tltle_node.getText())
    
    uls=soup.find_all('ul',class_='num_right')
    for ul in uls:
        for each in ul.find_all('li'):
            li_list=each.find_all('a')
            fund_info_dict={'fund_id':'',
                            'fund_name':'',
                            'fund_url':''}
            if len(li_list)>1:
                fund=li_list[0].text
                fund_id=re.findall(r'\d+',fund)[0]
                fund_url=li_list[0].attrs['href']
                fund_name=fund.decode('utf-8')[fund.find(ur') ') + 1:]encode('utf8')
                fund_info_dict['fund_id']=fund_id
                fund_info_dict['fund_name']=fund_name
                fund_info_dict['fund_url']=fund_url
                funds_text.append(fund_info_dict)
                return funds_text
    
class Handel_Url(Thread):
    All_Funds = []
    def __int__(self, queue):
        super(Handel_Url, self).__init__()
        self.queue=quue
        
