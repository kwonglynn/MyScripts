# -*- coding: utf-8 -*-
"""
Created on Thu May  2 10:28:59 2019

@author: Guanglin Kuang
"""

import scrapy
import json

class QuotesSpider(scrapy.Spider):
    name = "dict"

    def start_requests(self):
        url = 'https://www.dict.com/swedish-english/'
        word = getattr(self, 'word', None)
        if word is not None:
            url = url + word
        yield scrapy.Request(url, self.parse)

    def parse(self, response):
        entry_path = response.xpath(".//div[@id='entry-wrapper']/table[@class='entry']")
        word_path = entry_path.xpath(".//span[@class='lex_ful_entr l1']/text()")
        pron_path = entry_path.xpath(".//span[@class='lex_ful_pron']/text()")
        morf_path = entry_path.xpath(".//span[@class='lex_ful_morf']/text()")[0]
        
        trans_all = []
        trans = []
        phrase_all = []
        phrase = []
        for term in entry_path.xpath(".//td//span"):
            if term.attrib['class'] == 'lex_ful_labl':
                # Start append from the second translation.
                if trans != []:
                    trans_all.append(trans) 
                    trans = []
                                   
                trans.append(term.xpath("text()").get())
            elif term.attrib['class'] == 'lex_ful_prag':
                trans.append(term.xpath("text()").get())
            elif term.attrib['class'] == 'lex_ful_tran w l2':
                trans.append(term.xpath("text()").get())
            elif term.attrib['class'] == 'lex_ful_d':
                trans.append(term.xpath("text()").get())
                
            if term.attrib['class'] == 'lex_ful_coll2s w l1':
            # Start append from the second phrase.
                if phrase != []:
                    phrase_all.append(phrase)
                    phrase = []
                               
                phrase.append(term.xpath("text()").get())
            elif term.attrib['class'] == 'lex_ful_coll2t w l2':
                phrase.append(term.xpath("text()").get())
            
        # For the last translation
        trans_all.append(trans)                                          
        
        # For the last phrase
        phrase_all.append(phrase)
        
        yield {
            'word': word_path.get(),
            'pron': pron_path.get(),
            'morf': morf_path.get(),
            'trans': json.dumps(trans_all),
            'phrase': json.dumps(phrase_all)
        }
