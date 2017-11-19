/****************************************************************************
** Meta object code from reading C++ file 'EpipolarWindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../EpipolarWindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'EpipolarWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_EpipolarWindow_t {
    QByteArrayData data[7];
    char stringdata0[99];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_EpipolarWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_EpipolarWindow_t qt_meta_stringdata_EpipolarWindow = {
    {
QT_MOC_LITERAL(0, 0, 14), // "EpipolarWindow"
QT_MOC_LITERAL(1, 15, 20), // "epipolarLineOkSignal"
QT_MOC_LITERAL(2, 36, 0), // ""
QT_MOC_LITERAL(3, 37, 11), // "cv::Point3d"
QT_MOC_LITERAL(4, 49, 20), // "getLeftClickedPoints"
QT_MOC_LITERAL(5, 70, 21), // "getRightClickedPoints"
QT_MOC_LITERAL(6, 92, 6) // "accept"

    },
    "EpipolarWindow\0epipolarLineOkSignal\0"
    "\0cv::Point3d\0getLeftClickedPoints\0"
    "getRightClickedPoints\0accept"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_EpipolarWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    2,   34,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       4,    2,   39,    2, 0x0a /* Public */,
       5,    2,   44,    2, 0x0a /* Public */,
       6,    0,   49,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Int, 0x80000000 | 3,    2,    2,

 // slots: parameters
    QMetaType::Void, QMetaType::Double, QMetaType::Double,    2,    2,
    QMetaType::Void, QMetaType::Double, QMetaType::Double,    2,    2,
    QMetaType::Void,

       0        // eod
};

void EpipolarWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        EpipolarWindow *_t = static_cast<EpipolarWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->epipolarLineOkSignal((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< cv::Point3d(*)>(_a[2]))); break;
        case 1: _t->getLeftClickedPoints((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 2: _t->getRightClickedPoints((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 3: _t->accept(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (EpipolarWindow::*_t)(int , cv::Point3d );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&EpipolarWindow::epipolarLineOkSignal)) {
                *result = 0;
                return;
            }
        }
    }
}

const QMetaObject EpipolarWindow::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_EpipolarWindow.data,
      qt_meta_data_EpipolarWindow,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *EpipolarWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *EpipolarWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_EpipolarWindow.stringdata0))
        return static_cast<void*>(const_cast< EpipolarWindow*>(this));
    return QDialog::qt_metacast(_clname);
}

int EpipolarWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
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

// SIGNAL 0
void EpipolarWindow::epipolarLineOkSignal(int _t1, cv::Point3d _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
