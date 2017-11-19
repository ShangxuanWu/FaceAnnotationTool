/****************************************************************************
** Meta object code from reading C++ file 'VTKViewer.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../VTKViewer.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'VTKViewer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_VTKViewer_t {
    QByteArrayData data[16];
    char stringdata0[242];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_VTKViewer_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_VTKViewer_t qt_meta_stringdata_VTKViewer = {
    {
QT_MOC_LITERAL(0, 0, 9), // "VTKViewer"
QT_MOC_LITERAL(1, 10, 7), // "setHint"
QT_MOC_LITERAL(2, 18, 0), // ""
QT_MOC_LITERAL(3, 19, 8), // "point_id"
QT_MOC_LITERAL(4, 28, 9), // "setHint3D"
QT_MOC_LITERAL(5, 38, 10), // "cancelHint"
QT_MOC_LITERAL(6, 49, 12), // "cancelHint3D"
QT_MOC_LITERAL(7, 62, 17), // "meshPointerSignal"
QT_MOC_LITERAL(8, 80, 28), // "vtkSmartPointer<vtkPolyData>"
QT_MOC_LITERAL(9, 109, 17), // "to2DButtonPressed"
QT_MOC_LITERAL(10, 127, 29), // "correctThisPointButtonPressed"
QT_MOC_LITERAL(11, 157, 12), // "debugSetMesh"
QT_MOC_LITERAL(12, 170, 8), // "onCancel"
QT_MOC_LITERAL(13, 179, 26), // "enableFineAnnotationButton"
QT_MOC_LITERAL(14, 206, 16), // "toggleVisibility"
QT_MOC_LITERAL(15, 223, 18) // "to2DButtonPressed_"

    },
    "VTKViewer\0setHint\0\0point_id\0setHint3D\0"
    "cancelHint\0cancelHint3D\0meshPointerSignal\0"
    "vtkSmartPointer<vtkPolyData>\0"
    "to2DButtonPressed\0correctThisPointButtonPressed\0"
    "debugSetMesh\0onCancel\0enableFineAnnotationButton\0"
    "toggleVisibility\0to2DButtonPressed_"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_VTKViewer[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       7,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   74,    2, 0x06 /* Public */,
       4,    1,   77,    2, 0x06 /* Public */,
       5,    0,   80,    2, 0x06 /* Public */,
       6,    0,   81,    2, 0x06 /* Public */,
       7,    1,   82,    2, 0x06 /* Public */,
       9,    0,   85,    2, 0x06 /* Public */,
      10,    1,   86,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      11,    0,   89,    2, 0x0a /* Public */,
      12,    0,   90,    2, 0x0a /* Public */,
      13,    0,   91,    2, 0x0a /* Public */,
      14,    0,   92,    2, 0x0a /* Public */,
      15,    0,   93,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 8,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void VTKViewer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        VTKViewer *_t = static_cast<VTKViewer *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->setHint((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->setHint3D((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->cancelHint(); break;
        case 3: _t->cancelHint3D(); break;
        case 4: _t->meshPointerSignal((*reinterpret_cast< vtkSmartPointer<vtkPolyData>(*)>(_a[1]))); break;
        case 5: _t->to2DButtonPressed(); break;
        case 6: _t->correctThisPointButtonPressed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: _t->debugSetMesh(); break;
        case 8: _t->onCancel(); break;
        case 9: _t->enableFineAnnotationButton(); break;
        case 10: _t->toggleVisibility(); break;
        case 11: _t->to2DButtonPressed_(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (VTKViewer::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&VTKViewer::setHint)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (VTKViewer::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&VTKViewer::setHint3D)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (VTKViewer::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&VTKViewer::cancelHint)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (VTKViewer::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&VTKViewer::cancelHint3D)) {
                *result = 3;
                return;
            }
        }
        {
            typedef void (VTKViewer::*_t)(vtkSmartPointer<vtkPolyData> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&VTKViewer::meshPointerSignal)) {
                *result = 4;
                return;
            }
        }
        {
            typedef void (VTKViewer::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&VTKViewer::to2DButtonPressed)) {
                *result = 5;
                return;
            }
        }
        {
            typedef void (VTKViewer::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&VTKViewer::correctThisPointButtonPressed)) {
                *result = 6;
                return;
            }
        }
    }
}

const QMetaObject VTKViewer::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_VTKViewer.data,
      qt_meta_data_VTKViewer,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *VTKViewer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *VTKViewer::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_VTKViewer.stringdata0))
        return static_cast<void*>(const_cast< VTKViewer*>(this));
    return QWidget::qt_metacast(_clname);
}

int VTKViewer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 12)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void VTKViewer::setHint(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void VTKViewer::setHint3D(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void VTKViewer::cancelHint()
{
    QMetaObject::activate(this, &staticMetaObject, 2, Q_NULLPTR);
}

// SIGNAL 3
void VTKViewer::cancelHint3D()
{
    QMetaObject::activate(this, &staticMetaObject, 3, Q_NULLPTR);
}

// SIGNAL 4
void VTKViewer::meshPointerSignal(vtkSmartPointer<vtkPolyData> _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void VTKViewer::to2DButtonPressed()
{
    QMetaObject::activate(this, &staticMetaObject, 5, Q_NULLPTR);
}

// SIGNAL 6
void VTKViewer::correctThisPointButtonPressed(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
