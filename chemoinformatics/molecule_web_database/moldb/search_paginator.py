from math import ceil
from django.db.models import Q

class SearchPaginator(object):
    def __init__(self, model_class, ids_list, per_page, id_col_name="id"):
        try:
            self.__per_page = int(per_page)
        except (TypeError, ValueError):
            raise self.PageNotAnInteger("That get_page number is not an integer.")
        self.__model_class = model_class
        self.__object_ids_list = ids_list
        self.__objects_count = len(self.__object_ids_list)
        self.__pages_ids = self.__generate_pages_ids()
        self.__pages_count = len(self.__pages_ids)
        self.__id_col_name = id_col_name

    def get_page(self, number, return_ids=None):
        try:
            number = int(number)
        except (TypeError, ValueError):
            raise self.PageNotAnInteger("That get_page number is not an integer.")

        if number < 1:
            raise self.PageDoesNotExist(number)

        number -= 1

        try:
            ids = self.__pages_ids[number]
            if return_ids:
                return ids
            else:
                return self.Page(self.__get_objects(ids), number, self.__pages_count)
        except IndexError:
            raise self.PageDoesNotExist(number + 1)

    @property
    def object_ids(self):
        return self.__object_ids_list

    @property
    def objects_count(self):
        return self.__objects_count

    @property
    def pages_count(self):
        return self.__pages_count

    def __generate_pages_ids(self):
        pages = []
        for i in range(ceil(len(self.__object_ids_list) / self.__per_page)):
            start = i * self.__per_page
            end = start + self.__per_page
            try:
                pages.append(self.__object_ids_list[start:end])
            except IndexError:
                pages.append(self.__object_ids_list[start:])
        return pages

    def __get_objects(self, ids):
        return self.__model_class.objects.filter(Q((self.__id_col_name + "__in", ids)))

    class Page(object):
        def __init__(self, object_list, page_number, paginator_pages_count):
            self.__page_number = page_number + 1
            self.__paginator_pages_count = paginator_pages_count
            self.__objects_list = object_list
            self.__objects_count = len(object_list)

            if 1 < self.__page_number < paginator_pages_count:
                self.__has_next = True
                self.__has_previous = True
                self.__next_page_number = self.__page_number + 1
                self.__previous_page_number = self.__page_number - 1
            elif self.__page_number == 1 and paginator_pages_count == 1:
                self.__has_next = False
                self.__has_previous = False
                self.__next_page_number = None
                self.__previous_page_number = None
            elif self.__page_number == 1:
                self.__has_next = True
                self.__has_previous = False
                self.__next_page_number = self.__page_number + 1
                self.__previous_page_number = None
            elif self.__page_number == paginator_pages_count:
                self.__has_next = False
                self.__has_previous = True
                self.__next_page_number = None
                self.__previous_page_number = self.__page_number - 1

        def __getitem__(self, index):
            return self.__objects_list[index]

        @property
        def has_next(self):
            return self.__has_next

        @property
        def has_previous(self):
            return self.__has_previous

        @property
        def next_page_number(self):
            return self.__next_page_number

        @property
        def previous_page_number(self):
            return self.__previous_page_number

        @property
        def page_number(self):
            return self.__page_number

        @property
        def paginator_pages_count(self):
            return self.__paginator_pages_count

        @property
        def objects_count(self):
            return self.__objects_count

        @property
        def objects_list(self):
            return self.__objects_list

    class PageNotAnInteger(Exception):
        pass

    class PageDoesNotExist(Exception):
        def __init__(self, page_number):
            self.page_number = page_number

        def __str__(self):
            return "Page {} does not exist.".format(self.page_number)