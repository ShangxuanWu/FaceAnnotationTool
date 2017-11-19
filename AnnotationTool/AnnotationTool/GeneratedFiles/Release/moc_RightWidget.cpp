/****************************************************************************
** Meta object code from reading C++ file 'RightWidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../RightWidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'RightWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_RightWidget_t {
    QByteArrayData data[7];
    char stringdata0[64];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_RightWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_RightWidget_t qt_meta_stringdata_RightWidget = {
    {
QT_MOC_LITERAL(0, 0, 11), // "RightWidget"
QT_MOC_LITERAL(1, 12, 7), // "setHint"
QT_MOC_LITERAL(2, 20, 0), // ""
QT_MOC_LITERAL(3, 21, 8), // "point_id"
QT_MOC_LITERAL(4, 30, 9), // "setHint3D"
QT_MOC_LITERAL(5, 40, 10), // "cancelHint"
QT_MOC_LITERAL(6, 51, 12) // "cancelHint3D"

    },
    "RightWidget\0setHint\0\0point_id\0setHint3D\0"
    "cancelHint\0cancelHint3D"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_RightWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   34,    2, 0x0a /* Public */,
       4,    1,   37,    2, 0x0a /* Public */,
       5,    0,   40,    2, 0x0a /* Public */,
       6,    0,   41,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void RightWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        RightWidget *_t = static_cast<RightWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->setHint((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->setHint3D((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->cancelHint(); break;
        case 3: _t->cancelHint3D(); break;
        default: ;
        }
    }
}

const QMetaObject RightWidget::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_RightWidget.data,
      qt_meta_data_RightWidget,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *RightWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *RightWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_RightWidget.stringdata0))
        return static_cast<void*>(const_cast< RightWidget*>(this));
    return QWidget::qt_metacast(_clname);
}

int RightWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
